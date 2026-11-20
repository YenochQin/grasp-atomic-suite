# RMCDHF 的 MCP 数据需求与生成阶段预收缩

日期：2026-11-16。状态：**只读源码与数学依赖分析；未修改计算代码，未运行新的原子计算、
矩阵收缩或性能实验。** 源码基准：suite `caaf3a51ce12ec8e85e30622fc854c9e2a070683`；
上游参考 `63f2e1b29975337a136bcd95eb985ffef8507e89`。
范围是普通 `rmcdhf{,_mem,_mpi,_mem_mpi}`，以普通 MPI 版追踪主流程。
本篇续接 [MCP 存储分析](rangular_mpi_mcp_storage_2026-11-16.md)，重点是减少要保存的数学数据，
而非对原始字节做通用压缩。

## 结论

**可以把部分计算前移以减少数据，但不能将一般 EOL 所有 MCP 在生成时求和成一张固定势系数表。**

- 对通用、与径向函数和目标态无关的 MCP，可研究重复贡献归并、角向模板及索引映射。
  只要后续需要一般 EOL，就必须保留恢复逐 CSF 对算符作用的能力。
- 对冻结集合确定的 LBL 阶段，**所有参与轨道都冻结的积分贡献可以先算成固定矩阵 `H_fixed`**，
  替代这一部分展开的角系数列表。它在 EOL 混合系数更新后仍然有效；这是最值得优先量化的任务专用方案。
- 固定密度权重下可将角系数预收缩为每个积分一项；一般 EOL 的密度权重变化，不能一次收缩后丢弃原始映射。
- EOL 每次更新混合系数后，先做一次角系数收缩并缓存，能减少重复读取和计算；
  这是迭代期优化，单独使用不会减少永久 MCP 的必要信息量。

这里“数学等价收缩”与“逐位无损编码”不同：提前求和改变浮点运算顺序，
即使不引入物理近似，也不自动保证与旧程序逐位一致。

## 1. RMCDHF 实际消费的两种信息

| 使用阶段 | 需要 MCP 中什么 | 还依赖什么 | 输出 |
| --- | --- | --- | --- |
| MATRIX/SETHAM 装配矩阵 | block、稀疏位置、积分标签、每个 CSF 对的 T/V 系数 | 当前径向函数、网格、核势 | H 的各稀疏矩阵元，再对角化得 C |
| SETCOF 构造轨道方程系数 | 积分标签、CSF 对、T/V 系数 | 当前 C、目标态权重、广义占据 UCF | 较小的 DA/YA/XA 表 |
| YPOT/XPOT/DACON | 已收缩的 DA/YA/XA，不再直接读 MCP | 当前径向函数及其导数 | 网格上的势和非齐次项 |
| SETLAG | 经 SETCOF 等取得的势 | 轨道、冻结选择、占据 | 正交约束的乘子 |

依据：[SETHAM](../../src/appl/rmcdhf_mpi/setham.f90)、
[SETCOF](../../src/appl/rmcdhf_mpi/setcof.f90)、
[COFPOTMPI](../../src/appl/rmcdhf_mpi/cofpotmpi.f90)、
[YPOT](../../src/appl/rmcdhf_mpi/ypot.f90)、[XPOT](../../src/appl/rmcdhf_mpi/xpot.f90)、
[DACON](../../src/appl/rmcdhf_mpi/dacon.f90)、[SETLAGMPI](../../src/appl/rmcdhf_mpi/setlagmpi.f90)。

`.30` 提供的稀疏结构用于矩阵装配/对角化，也将系数的 `INDEX` 映射为 CSF 对。
SETHAM 丢弃读入的 `ICLMN`，但 SETCOF 用 `ICLMN` 和 `IROW(INDEX)` 得到两个 CSF 编号。
因此物理上需要的是 CSF 对身份，不是必须保持两个冗余索引字段。
四种程序需要的数学信息相同；内存版改变读取方式，未使 EOL 的依赖消失。

## 2. 用同一个角系数映射理解两条路径

令 p 表示 `(block,r,s)` 的单三角 CSF 对，q 表示某个一电子或二电子积分标签；
`A[p,q]` 为相应角系数。重复贡献如存在，可先理解成数学上的和。
`R[q;P,Q]` 是当前径向函数对应的积分。显式 MCP 部分给出：

```text
H[p] = H_analytic[p] + Σq A[p,q] R[q;P,Q]
```

`H_analytic` 包含由占据和 FCO/GCO 等直接计算、不需要显式 MCP 的对角贡献，
并非固定不变的矩阵。SETHAM 当前实际执行 `EMT(INDX) += TEGRAL*COEFF`。
RINTI/SLATER 读取当前 PF/QF；不能在 rangular 时用初始轨道算一次、以后一直沿用所有积分。
依据：[SETHAM](../../src/appl/rmcdhf_mpi/setham.f90)、
[RINTI](../../src/lib/core/one_electron_radial_integral.f90)、[SLATER](../../src/lib/core/slater_radial_integral.f90)。

EOL 的 CSF 密度权重由 DSUBRS 给出：

```text
D[b,r,s] = Σν WT[b,ν] C[b,r,ν] C[b,s,ν]
S[q] = Σb Σr≥s (2-δrs) D[b,r,s] A[b,r,s,q]
```

WT 在 NEWCO 中归一化。SETCOF 对每个标签先算上式的收缩，再乘 `0.5/UCF(J)`，
按目标轨道 J 在积分标签中出现的位置/次数分配到 DA、YA、XA，并合并不同 block 的同类项。
FCO/GCO 产生的隐式对角项另行按 D 的对角元素求和，不能漏掉。
依据：[DSUBRS](../../src/appl/rmcdhf_mpi/dsubrs.f90)、
[NEWCOMPI](../../src/appl/rmcdhf_mpi/newcompi.f90)、[SETCOF](../../src/appl/rmcdhf_mpi/setcof.f90)。

这说明 A 有两种用途：给定径向积分计算 H；给定密度权重计算 S。
保存 S 一般不能恢复 A，也不能计算下一轮所需的 H 或不同 D 对应的 S。
把 MCP 改称“逐 CSF 对的跃迁密度矩阵”也不会自动减少其信息量。

例如同一积分 R 的两个对角角系数为 1 和 −1：

```text
H = R × diag(1,-1)
```

跨 CSF 对提前相加得到 0，会丢失本来的两个本征值 ±R。
只有先知道所需的固定 D，才能决定这两项应如何加权；也不能把不同径向积分的角系数直接相加。

## 3. 一般 EOL 为何不能一次预收缩到底

[SCFMPI](../../src/appl/rmcdhf_mpi/scfmpi.f90) 在开始时及每轮轨道更新后调用
`MATRIXMPI → NEWCOMPI`。矩阵对角化更新 EVEC，NEWCO 重算 UCF；
这些量随后进入 SETCOF。即使优化态的**权重 WT 固定，C 仍在变化，所以 D 不固定**。
冻结径向轨道也不意味着冻结 ASF 的 CSF 混合系数。

只保留初始态的 S 会将后续计算约束到旧密度，改变一般 EOL 的算法。
只保留目标 ASF 的初始投影同样不能保留对角化过程中的全部耦合。
选定一个 J/宇称块里的少数目标根，并不意味着只需保存这些根初猜中较大的 CSFs。

有严格条件的例外：固定 C 的轨道优化、正确实现的固定权重 EAL，或完整子空间等权迹平均。
即便轨道目标函数的 D 固定，若还要求一般 H 的本征值/本征向量，仍需另保留构造 H 的能力。

## 4. 最值得量化的缩减：冻结积分收缩为 H_fixed

冻结集合 F 必须指**本阶段径向函数真正不再变化的轨道集合**，不是“禁止从它激发”的轨道集合。
径向冻结和占据固定是两件事：冻结轨道仍可在不同 CSF 中有不同占据和耦合。

将积分标签分为两类：

| 标签 | 全冻结条件 | 处理 |
| --- | --- | --- |
| I(a,b) | a、b 均在 F | 径向积分可预计算为常数 |
| Rᵏ(a,b,c,d) | a、b、c、d 均在 F | 径向积分可预计算为常数 |
| 其余 | 至少一个轨道会变化 | 仍保留角系数与可变积分的对应关系 |

对全冻结集合 q∈Fint，定义：

```text
H_fixed[p] = Σq∈Fint A[p,q] R_fixed[q]
H[p;P,Q] = H_fixed[p] + H_analytic[p;P,Q]
           + Σq∉Fint A[p,q] R[q;P,Q]
```

这样原来同一 CSF 对上许多全冻结积分的角系数条目，变为一个固定矩阵元素。
H_fixed 不依赖 C，因此 **EOL 继续更新混合系数并不使它失效**。
它仍影响新本征向量，必须加入每次矩阵装配；不能把整个 H_fixed 简化成一个全局常数。
只有另外证明某部分在所选 CSF 空间严格等于 `E_core × I`，才能单独保存该标量。

对可变轨道 J，全冻结积分不含 J，故其对 J 的径向变分为零；SETCOF 本来就按标签是否包含 J 跳过。
这给出在活性轨道势构造时去掉这些展开条目的依据。
隐式 FCO/GCO/占据项可保持原来计算，也可另外分类，但不能既加入 H_fixed 又重复计算。

### 对 SETLAG 的精确限制

当前四种普通程序都会先构造一对轨道两侧的 SETCOF/势，即使其中一侧固定。
但一固定一可变时，最终乘子公式只使用**可变侧**的势：
普通 MPI 的 `LFIX(M)` 分支用 J 侧，`LFIX(J)` 分支用 M 侧。
两轨道均固定的对不建立乘子；两轨道均可变时需要两侧势。

所以不能仅凭“冻结侧也调用 SETCOF”就判定全冻结积分必须永久保留；
固定侧的部分工作目前是冗余的。但是直接删旧文件中的条目也不是现有程序已支持的路径。
将来必须调整或兼容这些访问，且先核对固定/可变分支、原占据因子和 MPI 归约。
依据：[普通 SETLAGMPI](../../src/appl/rmcdhf_mpi/setlagmpi.f90)、
[内存 SETLAGMPI](../../src/appl/rmcdhf_mem_mpi/setlagmpi.f90)、
[串行 SETLAG](../../src/appl/rmcdhf/setlag.f90)、
[内存串行 SETLAG](../../src/appl/rmcdhf_mem/setlag.f90)。

### 收益怎样估算

设被替换的全冻结角系数共有 Nf 条，它们落在 Mf 个不同稀疏位置上。
原贡献主体约 `16 Nf` 字节；若新固定矩阵保存 `(int32 INDEX,float64 value)`，
主体约 `12 Mf` 字节，另外有 block/目录开销，并保留所需 `.30` 结构。
若复用完整结构存密集位置数组，则是 `8 NNONZ` 字节；二者应按实际分布选择。
**收益取决于 Nf、Mf，而非冻结轨道数量本身**。没有实际冻结列表，本轮不报缩减百分比。

对逐层优化、以前各层径向轨道全部冻结的阶段，这种任务专用收缩值得优先统计。
若之后解冻轨道、换 rwfn、网格、核势或算符设置，相关 H_fixed 必须失效并重建。
轨道相位/归一化及积分精度约定也须固定；不能只凭 `LFIX` 标记断言函数不变。
标准 EOL 的 [ORTHY](../../src/appl/rmcdhf_mpi/orthy.f90) 按固定/可变轨道分段处理，
而 `ORTHST=false` 时调用的 [ORTHSC](../../src/lib/core/radial_orbital_orthogonalization.f90) 没有 LFIX 过滤，
可能通过正交化改变标为固定的函数；这种路径必须另外检查。
重新求和可能引入舍入差异：需验证矩阵和完整 SCF 行为，不能称为已通过逐位一致验证。

### 应放在哪一步

当前 [rangular_mpi](../../src/appl/rangular_mpi/genmcpmpi.f90) 不读取 rwfn 和 RMCDHF 的冻结选择；
这些信息在 [GETSCDMPI](../../src/appl/rmcdhf_mpi/getscdmpi.f90) 等 RMCDHF 输入阶段才确定。
因此有两种不同设计：

1. 保留通用 rangular，在 RMCDHF 启动时流式构造 H_fixed，之后使用缩减存储。
   侵入生成端较少，但原始完整 MCP 仍曾落盘，不能消除首次生成峰值。
2. 让生成/预处理流程显式接收冻结列表、固定 rwfn、网格和模型参数，逐 CSF 对直接写 H_fixed 与可变项。
   可减少首次最终输出，但产物成为“本次冻结阶段专用数据”，不再是任意径向函数都适用的纯角系数库。

无论采用哪一种，均需同步适配消费端。当前未实施任何一种。

## 5. 与径向函数无关的进一步计算

### 同一 CSF 对、同一积分标签的归并

如果内部角动量求和给出多条相同 `(block,r,s,q)`，数学上可在生成阶段加成一个系数。
不能跨不同 CSF 对、不同 k 或不同积分归并；相同数值也不代表相同物理贡献。
前一轮 38,461,197 条贡献样本中已经没有这种重复，故对该样本没有收益；
其他模型须重新统计。浮点提前求和亦可能改变逐位结果。

### 角向模板与映射

角向等价的 CSF 对可以共享一份系数模板，再保存真实 CSF 对及径向轨道标签的映射。
后续 H 装配和 S 收缩使用模板展开的算符作用，不将多个 CSFs 合成一个混合系数。
这是仍可支持一般 EOL 的结构性缩减方向。

等价判断必须保留 κ、占据、子壳层状态、完整耦合链、轨道顺序、bra/ket 的共享轨道身份和相位。
不能只用相同 l/j 或省去主量子数就直接认定可复用。
本仓库已有串行 [rangular_reuse](../../src/appl/rangular/rangular_reuse.f90) 的保守键与映射，
但目前命中后仍展开写原 MCP，**只复用计算，不缩减输出**；MPI 未接入。
[其文档](../implemented/rangular_angular_reuse.md)及[历史测量](rangular_reuse_benchmark_2026-11-14.md)
还显示收益强烈依赖数据、有限缓存和遍历顺序，不能把命中率当作新格式压缩率。

CSFG 文献提供了“保存生成元的角系数，通过映射恢复完整 CSF 算符作用”的已有路线，
其 MCDHF 仍求完整展开空间的混合系数。可参考其信息组织思想，但不能据此宣布任意截断/ML 选择后的
CSF 列表都能直接套用完整族规则；缺失成员必须准确表达。
本轮未启用或迁移 CSFG。[原始论文，第 4.4–4.5 节](https://arxiv.org/html/2410.11297v1#S4.SS4)。

### 闭壳层与部分冻结积分

当前程序已经把一些简单对角项排除在 MCP 外：生成端调用 `RKCO_GG(...,INCOR=0,...)`、
排除 F0，消费端用占据及 FCO/GCO 计算相应部分。
不能再次把“所有闭壳层项都改成公式”当成未利用的新容量收益。
依据：[MCPMPI](../../src/appl/rangular_mpi/mcpmpi_gg.f90)、
[RKCO_GG](../../src/lib/spin_angular/two_electron_angular_coefficients.f90)、
[FCO](../../src/appl/rmcdhf_mpi/fco.f90)、[GCO](../../src/appl/rmcdhf_mpi/gco.f90)。

只冻结二电子积分中的一对轨道，可预先计算相应径向核，或在额外证明完全闭壳层条件下构造有效算符；
但剩余可变轨道仍决定积分，不能直接变成一个静态矩阵数值。
保存每个组合的径向网格函数还可能比保存几个角系数更大，需另做因子化/存储成本分析。

## 6. 固定密度的特例与当前 EAL 边界

若 D 真正固定且只要求该固定密度的轨道目标函数，则可以预先保存每个积分的 S[q]，
加上隐式对角贡献的固定加权系数；不再保存每个 CSF 对的展开表。
EAL 的理论形式为 `D_rs=δ_rs*w_r`，这是比一般 EOL 强得多的约束。
固定目标态权重 WT、固定轨道、固定 C 三者不能混淆。

当前源码不支持把它立即当作可靠现成方案：`GETSCDMPI` 在 DIAG 为真时选择 EAL，
其他普通路径选择 EOL；SCF 的 EAL 分支存在尚未验证的状态初始化、block 权重寻址和能量判据问题。
特别是 DSUBRS 的 EAL 分支直接取 `WT(I)`，未使用 block 偏移；GETALD 未走 GETOLD 的 NEVBLK
等初始化，而 SCF 在 EOL 条件分配之前使用这些 block 数组。
这些是静态风险，不是本轮运行复现的失败，不能据此推销切到 EAL 作为压缩办法。
依据：[GETSCDMPI](../../src/appl/rmcdhf_mpi/getscdmpi.f90)、
[GETALDMPI](../../src/appl/rmcdhf_mpi/getaldmpi.f90)、
[GETOLDMPI](../../src/appl/rmcdhf_mpi/getoldmpi.f90)、
[SCFMPI](../../src/appl/rmcdhf_mpi/scfmpi.f90)、[DSUBRS](../../src/appl/rmcdhf_mpi/dsubrs.f90)。

## 7. 每次混合系数更新后的收缩：省读取但不替代通用 MCP

SETCOF 当前对不同轨道 J 重复读取，同一标签的 `Σ D*A` 在计算不同 J 时重复出现。
在一次 MATRIX/NEWCO 更新之后，可先遍历各标签计算 S[q]，再由标签派生各 J 的 DA/YA/XA。
在 C、WT、UCF 不变期间，这些角向收缩结果可以复用；下一次更新后必须失效。
隐式对角项也可按同一密度版本缓存，MPI 可保留局部收缩或对小表归约，但不能重复计数。

**可缓存的是收缩系数，不是无条件缓存整个 YP/XP 网格势。** 后者仍依赖正在逐个改变的径向轨道。
这类改法可能减少 SCF 期间的 MCP 扫描与 DSUBRS 计算，但下一轮仍需原始 A 或等价的模板/因子化表示。
因此应与永久数据缩减分别测量，不能混报容量收益。

## 8. 下一步只读分析应测什么

在决定是否实现前，针对实际 LBL 任务给定冻结列表与目标 block，统计：

1. 全冻结标签的贡献数 Nf、覆盖位置数 Mf，各 block/k 的占比；估计 H_fixed 与剩余可变项总量。
2. 不含目标根且无其他用途的整个 J/宇称 block 的量；这些块可在任务专用数据中排除，
   但不能删掉被优化块内部未直接选作目标态的 CSFs。
3. 同一 `(CSF对,积分)` 的重复率，以及完整角向模板的数量、模板长度和映射成本。
4. 每轮 NEWCO 后 SETCOF 扫描次数与各标签复用次数，区分容量瓶颈和读取/计算瓶颈。

本轮没有给任意冻结掩码代入真实样本，也没有构造 H_fixed；故没有声称缩减倍数或数值验收。
推荐先量化“**全冻结积分 → H_fixed；含可变轨道积分 → 保留角系数**”，
再决定是否值得引入模板存储和迭代期收缩缓存。通用全变轨道 EOL 则应优先研究模板/编码，
不能依赖一次密度预收缩取代所需的算符信息。
