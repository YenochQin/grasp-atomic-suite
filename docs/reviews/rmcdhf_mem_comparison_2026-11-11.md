# RMCDHF 内存版与普通版源码对比

核查日期：2026-11-11。suite 源码基准：`db3950563fb9790fb5dccf609ccd51d766949496`。
本地 `grasp` 参考树基准：`ad3fb76636cab47f08b2accc8e6fb9d1650dc3b2`。

范围：普通 `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi`。
状态：已完成源码静态对比；未修改计算实现，未运行新的原子计算或性能基准。
不涉及 CSFG 或已退役的 orbopt、strict、diag 变体。

**`mem` 的主要作用是将 MCP 角系数和索引缓存到内存，但当前两类实现还存在
拉格朗日乘子筛选、MPI 乘子公式及径向边界处理差异。不能将它们视为只改变
存储位置、保证数值行为相同的两个执行方式。**

## 程序共同部分

串行两版的 `scf.f90` 逐字节相同，MPI 两版的 `scfmpi.f90` 逐字节相同。
相应的 `MANEIG`、`NEWCO`、`SOLVE`、`CONSIS`、阻尼和正交化例程也相同。
这表示对应版本沿用相同的外层 SCF 控制、对角化实现和停止判据，
不表示它们调用的全部下层例程相同。

按同名 `.f90` 文件比较，串行两版有 110 个共同文件，其中 94 个相同、16 个不同；
MPI 两版有 113 个共同文件，其中 97 个相同、16 个不同。每个 `mem` 目录另有
17 个文件，包含缓存模块、加载/读取例程、内存查询及接口文件。
文件统计包含接口文件，不能直接解释为算法差异数量。

依据：[串行 SCF](../../src/appl/rmcdhf/scf.f90)、
[内存版串行 SCF](../../src/appl/rmcdhf_mem/scf.f90)、
[MPI SCF](../../src/appl/rmcdhf_mpi/scfmpi.f90)、
[内存版 MPI SCF](../../src/appl/rmcdhf_mem_mpi/scfmpi.f90)。

## MCP 存储与读取

| 比较项 | 不带 mem | 带 mem |
| --- | --- | --- |
| 初始化 | 打开和检查 MCP 文件，供后续读取。 | 扫描记录以确定大小，分配缓存，再回读加载并关闭相应文件。 |
| 矩阵和轨道方程系数构造 | `MATRIX`、`SETHAM`、`SETCOF` 反复定位、读取 MCP。 | 使用 `read*_mem` 从缓存取出标签、系数及索引。 |
| 常驻数据 | 不额外保存完整 MCP 缓存。 | 保存全部所需 MCP 数据；MPI 版保存各进程对应的 MCP 数据。 |
| 仍需计算的量 | 根据当前轨道计算径向积分并重建 Hamiltonian。 | 同样计算径向积分并重建 Hamiltonian。 |

缓存内容包括：`mcp.30` 的稀疏矩阵结构 `IENDC`、`IROW`；`mcp.31` 的一电子
MCP 贡献；`mcp.(32+k)` 的二电子 MCP 贡献，以及各文件的 `LAB`、`NCONTR`、
`ICLMN`、`INDX`、`COEFF`。缓存对象是角系数和结构信息，
随轨道变化的径向积分与 Hamiltonian 数值没有因此固定。

`read2_mem`、`read3_mem` 仍将缓存内容复制到工作数组，`SETHAM` 等仍有工作区
分配和扩容；不能将内存版描述为消除了全部数据搬运或分配开销。

MPI 入口按 `'mcp' // idstring` 读取各 rank 的文件，并检查文件中记录的
`MYID/NPROCS` 与运行配置一致。它没有让每个 rank 加载一份完整的全局角系数；
也不能由此假定所有数组都按 rank 数等比例缩小。

依据：[SETMCP](../../src/appl/rmcdhf_mem/setmcp.f90)、
[SETMCP_MEM](../../src/appl/rmcdhf_mem/setmcp_mem.f90)、
[缓存数组](../../src/appl/rmcdhf_mem/rmcdhf_mem_C.f90)、
[READ2_MEM](../../src/appl/rmcdhf_mem/read2_mem.f90)、
[READ3_MEM](../../src/appl/rmcdhf_mem/read3_mem.f90)、
[MPI 入口](../../src/appl/rmcdhf_mem_mpi/rscfmpivu.f90)。

## 拉格朗日乘子轨道对筛选

串行和 MPI 的 `mem` 与非 `mem` 对比都存在此差异。
两类版本均只考虑同一 kappa、满足 `IORDER/NCORE` 遍历范围的候选轨道对，
并要求至少一个轨道参与优化。满占据判断为：

```fortran
EPS = ACCY*0.01D0
FULLI = ABS(UCF(LI)-DBLE(NKJ(LI)+1)) < EPS
```

不带 `mem` 的版本排除“两个轨道都优化且都满占据”的轨道对；
带 `mem` 的版本排除所有“两个轨道都满占据”的轨道对。

| 候选轨道对状态 | 不带 mem | 带 mem |
| --- | --- | --- |
| 均满占据，一个固定、一个优化 | 建立乘子 | 不建立乘子 |
| 均满占据，两个都优化 | 不建立乘子 | 不建立乘子 |
| 两个都固定 | 不建立乘子 | 不建立乘子 |
| 至少一个优化，且不同时满占据 | 建立乘子 | 建立乘子 |

`NEC/IECC` 在首次调用时建立，因此这里的满占据筛选使用首次建表时的 `UCF`。
此外，两类版本的嵌套遍历顺序不同，所列轨道对的顺序也可能不同。
筛选条件差异会改变实际建立的乘子集合，不能归因于缓存或浮点舍入。

依据：[普通版 SETLAG](../../src/appl/rmcdhf/setlag.f90)、
[内存版 SETLAG](../../src/appl/rmcdhf_mem/setlag.f90)、
[普通版 SETLAGMPI](../../src/appl/rmcdhf_mpi/setlagmpi.f90)、
[内存版 SETLAGMPI](../../src/appl/rmcdhf_mem_mpi/setlagmpi.f90)。

## MPI 拉格朗日乘子公式选择

对于已入选且两个轨道都参与优化的轨道对，内存 MPI 版计算：

```fortran
QDIF = ABS((UCFJ - UCFM)/MAX(UCFJ,UCFM))
```

变量 `P001` 的实际定义是 `1.0D-01`，即 **0.1**。
内存 MPI 版在 `QDIF > 0.1` 时，采用含势之差、分母为
`1/UCFJ - 1/UCFM` 的表达式；其余情况采用含势之和、分母为
`1/UCFJ + 1/UCFM` 的表达式，并包含对应的一电子积分修正。
普通 MPI 版在此分支直接采用后一种表达式，没有 `QDIF` 切换。

| 程序 | 两个可变轨道的公式选择 |
| --- | --- |
| `rmcdhf` | 根据 `QDIF > 0.1` 切换两种表达式。 |
| `rmcdhf_mem` | 根据 `QDIF > 0.1` 切换两种表达式。 |
| `rmcdhf_mpi` | 始终采用含势之和的表达式。 |
| `rmcdhf_mem_mpi` | 根据 `QDIF > 0.1` 切换两种表达式。 |

因此，此项是 **MPI 两版之间**的实质公式选择差异；不能把它写成串行两版
也存在同样的公式分支差异，也不能把 MPI 两版的全部迭代差异都归为归约舍入。
静态对比本身不裁定哪种处理在某个体系中更合适。

依据：[普通 MPI 公式](../../src/appl/rmcdhf_mpi/setlagmpi.f90)、
[内存 MPI 公式](../../src/appl/rmcdhf_mem_mpi/setlagmpi.f90)，
以及上一节的两个串行 `SETLAG`。

## 径向积分外边界

串行和 MPI 的 `mem` 与非 `mem` 对比都存在此差异。
不带 `mem` 的 `IN` 接受条件为：

```fortran
IF (TTHIS + TLAST <= EPS .OR. J .EQ. NNNP-1) THEN
   MTP = J
```

带 `mem` 的接受条件为：

```fortran
IF (TTHIS + TLAST <= EPS) THEN
   MTP = J
```

两版在推进时均先检查 `J >= N`，达到运行网格限制便停止。
在 `N = NNNP` 的情况下，如果衰减条件尚未满足，普通版仍可在 `NNNP-1`
处接受外边界；内存版继续推进后会触发网格限制报错。
若尾部提前满足衰减条件，这项分支差异不会被触发。

`NNNP` 是编译期容量，`N` 是运行时网格点数，二者不能混为一谈。
当 `N < NNNP` 时，不能假定普通版总能靠 `NNNP-1` 条件继续运行。
普通版没有报错也不足以证明轨道尾部已充分衰减。

依据：[普通串行 IN](../../src/appl/rmcdhf/in.f90)、
[内存串行 IN](../../src/appl/rmcdhf_mem/in.f90)、
[普通 MPI IN](../../src/appl/rmcdhf_mpi/in.f90)、
[内存 MPI IN](../../src/appl/rmcdhf_mem_mpi/in.f90)。

## 资源开销和异常处理

从实现可预期，反复扫描 MCP 且文件读取占比较高的计算可能受益于内存版；
若主要时间在径向积分、对角化或 MPI 通信中，缓存收益可能有限。
这是性能预期，尚无本次实测速度或峰值内存数据。

在默认 4 字节整数、8 字节实数的布局下，每条缓存贡献的
`ICLMN + INDX + COEFF` 约占 16 字节；每个标签及贡献数对另占约 8 字节。
总内存还包括稀疏结构、block 偏移、临时工作区、矩阵和本征向量，
不能只用角系数文件大小作为作业内存申请值。

`SYSTEMMem/SYSTEMFreeMem` 读取 Linux `/proc/meminfo` 并打印信息，
缓存加载代码未实现内存不足时自动回退到文件读取。
其显示的系统内存信息也不等同于调度器授予该作业的内存限额。

异常处理还有独立差异：普通版多处使用 `ERROR STOP`，内存版对应位置使用
`STOP`；内存版 `SETMCP` 的累积错误分支关闭文件后，缺少普通版对应的立即终止
语句。异常输入下的退出和后续行为需要单独验证，不能由正常输入路径推断。

依据：[缓存分配](../../src/appl/rmcdhf_mem/setmcp_mem.f90)、
[系统总内存查询](../../src/appl/rmcdhf_mem/systemmem.f90)、
[系统空闲内存查询](../../src/appl/rmcdhf_mem/systemfreemem.f90)、
[普通 SETMCP](../../src/appl/rmcdhf/setmcp.f90)、
[内存 SETMCP](../../src/appl/rmcdhf_mem/setmcp.f90)。

## 验证范围与后续对照

本次对比覆盖四个目录的同名 Fortran 文件及内存版新增文件。
上述三个数值差异涉及的四份 `SETLAG/SETLAGMPI` 和四份 `IN`，
分别与所列本地 `grasp` 参考版本中的对应文件逐字节一致；
因此这些差异在本地参考树中已有，不能认定为 suite 本次记录所引入。
此结论限定于所列版本和文件，不代表核查了远端最新上游。

源码证据确认了条件与公式不同，尚未量化它们对具体算例的迭代次数、
能量、轨道或其他可观测量的影响，也没有建立物理精度或活性空间收敛结论。

后续对照应固定 CSF、初猜、核模型、网格、优化轨道、权重、阻尼和编译设置；
MPI 两版还应使用相同 rank 数及对应的相同 MCP 输入。
先比较首次乘子轨道对、`UCF/QDIF/ECV`、径向 `MTP` 与各轮轨道和能量，
确认数值分支是否触发，再分别统计初始化、系数构造、对角化及总耗时和内存。
若要隔离纯缓存收益，应另行设计计算逻辑受控一致的对照，避免将不同迭代路径
造成的耗时变化误判为存储优化收益。此受控对照尚未实施。

相关说明：[普通串行 RMCDHF 流程](../implemented/rmcdhf_flow_and_state.md)、
[公共参数](../implemented/common_parameters.md)、
[单节点 MPI 诊断计划](../plans/rmcdhf_mpi_single_node_diagnostics_2026-11-09.md)。
