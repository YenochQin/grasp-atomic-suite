# GRASP 径向网格参数调整与验证规程

适用对象：使用 `scripts/patch_grasp_grid.py` 修改完整原版 GRASP，保留多份
独立源码目录，按原版 README 的 CMake 流程手动编译。另说明
`grasp-atomic-suite` 集中参数模式的差别。

编写日期：2026-10-01。脚本接口依据 suite 提交 `b60cbde`；原始默认值依据
本工作区参考 GRASP 提交 `fc81312614e3c1aff18090723a7d4e79260a49f1`。
服务器版本不同时，应先查看其源码和运行输出，不假定默认值完全相同。

本文给出可以执行的**参数试验方案**，不是已经完成的物理收敛报告。
原版手册提供扩展网格实例；具体体系是否合适，由本文规定的计算数据决定。

## 1. 必须先区分的四个问题

1. **容量**：`NNNP` 决定数组可以容纳多少点，`NNN1=NNNP+10`。
2. **网格分辨率**：主要由 `H` 和 `RNT` 决定。
3. **覆盖范围**：由 `N`、`H`、`RNT`、`HP` 共同决定。
4. **迭代与数值阈值**：`ACCY` 被数值程序及轨道自洽判断使用；不是物理量误差的保证。

因此，“容量改为 2990”与“网格已经足够精细”是两项不同结论。
脚本开头预填的 `nnnp=2990, n=2990, h=0.05` 是接口示例，尚未针对生产体系验证。
同样，suite 当前全部有限核默认 `N=NNNP=2990` 是统一配置策略，不能替代收敛试验。

**本方案的源码基底与原版基线均为未经修改的 GRASP，`NNNP=590、NNN1=600`。**
当前已调整到 2990 的源码或 suite 不能充当原版基线。保留一份原版 B0，
其他原版测试副本都从同一份干净原版复制；2990 仅是后续扩容组的候选容量。

## 2. 参数原值、候选值、原因和检查项

下面的“原值”是参考源码中的典型默认值。原版的不同程序存在历史差异；
脚本会同步所覆盖入口，但以运行时实际值为最终依据。

| 文件配置键 | 参考原值 | 本文候选值 | 调整原因 | 必须检查 |
| --- | --- | --- | --- | --- |
| `nnnp` | 590 | 有限核扫描 2990；不足时再试 3990；点核扫描 990 | 为较小 H 所需的实际点数提供容量；同一收敛系列固定容量 | 所有库、串行/MPI 程序全量重编译；`N<=NNNP`；内存和时间 |
| 派生 `NNN1` | 600 | 3000、4000 或 1000，脚本自动设置 | `RADGRD` 还生成 N 之后的 10 个辅助点 | 不独立修改；不能只改主模块而遗漏重复声明 |
| `n` | 有限核轨道程序通常为 `NNNP`，初始容量为 590 | 分辨率试验：590、1179、1965、2946 | 配合 H 调整，使外端范围大致相同 | `R(N)`、轨道 MTP、外层尾部及物理量 |
| `h` | 0.05 | 0.025、0.015、0.010 | 在相近范围内增加采样密度 | 能级间隔、轨道形状、性质积分；不能固定 N 直接减半 H |
| `rnt_scale` | `2e-6` | `2e-6 → 1e-6 → 5e-7` | 改变近核尺度，检验内层及核附近积分敏感性 | 同时补偿 N，避免把近核与截断误差混在一起 |
| `hp` | 0 | 第一轮保持 0 | 保持指数网格，便于控制变量 | 非零时须重算范围，不能套用指数网格公式 |
| `accy` | 默认公式 `H**6` | 网格扫描先固定 `1e-10`；随后 `1e-11、1e-12` | 防止 H 改变时同时改变求解阈值；另测阈值收敛 | 是否严格收敛、迭代次数、节点判断及最终物理量 |
| `point_n` | `MIN(220,NNNP)` | 点核单独试验：220、439、877 | 与点核 H 配套，保持覆盖范围 | 不与有限核模型混合比较 |
| `point_h` | 0.0625 | 0.03125、0.015625 | 点核分辨率加密 | 单电子基准及点核轨道积分 |
| `point_rnt_scale` | `exp(-65/16)`，约 0.01720595 | 第一轮 `None`，保留原值 | 与有限核尺度的差别很大，不宜直接套用有限核设置 | 明确当前核模型和实际 RNT |

`rnt_scale` 是 **RNT 的系数**：`RNT=rnt_scale/Z`，不是直接输入的实际 RNT。
`None` 表示脚本不修改对应数值；0 是数值本身，两者不能互换。
原版 `accy` 显式设置必须为正数；suite 可用 `accy=0` 恢复自适应 `H**6`。

`NNNW=127` 和 `KEYORB=215` 不属于本次网格调整，脚本不修改它们。
`KEYORB` 有整数打包约束，不能因为增加径向点数而跟着修改。
节点阈值 `THRESH` 也不是脚本参数：原版默认 0.05，suite 对应
`NODE_THRESHOLD=0.05`。网格试验先保持它不变，避免更换判据掩盖异常节点。
若单独诊断，可比较 0.025、0.05、0.10，并检查真实轨道形状；阈值变大只会
让较小振荡被忽略，不能据此判定轨道正确。该机制见
[COUNT 源码](../src/lib/lib9290/count.f90)。

## 3. 用实际公式解释为什么要这样改

### 3.1 HP=0：指数网格

根据 [RADGRD](../src/lib/lib9290/radgrd.f90)，采用原子单位：

```text
R(1) = 0
R(i) = RNT * [exp((i-1)*H) - 1]
RNT = rnt_scale/Z
Rmax = R(N)
第一个非零节点 R(2) = RNT * [exp(H)-1]
```

所以 RNT 本身不是这个实现的第一个非零节点。
在固定 H/RNT 下增加 N，原有节点位置不变，只延长网格。
在固定 N/RNT 下减小 H，网格变细但外端也会显著缩短。

例如 Z=92、`rnt_scale=2e-6`：

| N | H | `(N-1)H` | 计算得到的 R(N)，a0 | 意义 |
| --- | --- | --- | --- | --- |
| 590 | 0.05 | 29.450 | `1.34034e5` | 参考默认指数范围 |
| 1179 | 0.025 | 29.450 | `1.34034e5` | 同范围，加密 H |
| 1965 | 0.015 | 29.460 | `1.35381e5` | 约同范围，更细 H |
| 2946 | 0.010 | 29.450 | `1.34034e5` | 同范围，更细 H |
| 2990 | 0.05 | 149.450 | `1.74804e57` | 只增加 N 会生成极远的外端，原有区域没有加密 |

这些是按公式计算的**网格外端**，不是电子云实际半径。
轨道可能在较小 MTP 后就被截断；很大的 R(N) 既不证明精度高，也不保证实际
计算一定失败。应查看有效轨道区间、尾部和实际耗时。

给定目标 Rmax，可估算需要的 N：

```text
N = 1 + ceil(log(1 + Rmax/RNT)/H)
```

粗细网格保持相近范围时，可先按下面的关系选择整数 N：

```text
(N_new-1)*H_new ≈ (N_old-1)*H_old
```

由此得到本文的 590/1179/1965/2946 组合。若实际研究需要更大外端，应单独
做范围试验，不应把所有体系强制限制在这组默认范围内。

### 3.2 HP>0：另一个网格族

同一源码在 HP 非零时求解：

```text
log(1 + R/RNT) + (H/HP)*R = (i-1)*H
```

足够远处渐近为线性网格，相邻点间距趋近 HP。HP 的单位为 a0。
改变 HP 会同时改变远端布局；本文第一轮保持 HP=0。
若确需非零 HP，必须重新设计 N 和范围，并完整重复分辨率、尾部与性质验证。

### 3.3 H 与 ACCY 的耦合

参考原版 `src/appl/rmcdhf90/getscd.f90` 的默认公式为 `ACCY=H**6`；
suite 将同一公式集中在 [radial_grid_defaults](../src/lib/libmod/radial_grid_defaults_M.f90)。
例如：

| H | 自适应 ACCY |
| --- | --- |
| 0.05 | `1.5625e-8` |
| 0.025 | `2.44140625e-10` |
| 0.015 | `1.1390625e-11` |
| 0.010 | `1e-12` |

H 减半，默认阈值会缩小 64 倍。结果改善可能同时来自网格加密和更严格的
迭代。本文首先固定 ACCY，再单独验证 ACCY；若固定的 `1e-10` 未满足误差预算，
用通过验证的更严格固定值重跑相关网格，不能沿用被求解误差污染的结论。
ACCY 在 [SCF](../src/appl/rmcdhf90/scf.f90) 中参与 SCNSTY 判断，也在
数值求解和节点判断中使用，不能把它直接当成能量误差或所有性质的相对误差。

## 4. 手册的 U I 示例应该怎样理解

[GRASP Manual for Users，§13.5，PDF 印刷页 343、346–348](https://www.diva-portal.org/smash/get/diva2:1771324/FULLTEXT01.pdf#page=343)
使用 U I 的 `5f³6d7s²` 构型，Z=92、A=238，扩展容量到
`NNNP=1990、NNN1=2000`，实际输入 `N=1990、H=0.015、HP=0、RNT=2.17e-8`。
这说明扩容与减小 H 要配合进行；不是所有重元素都应复制同一数值。

通过脚本近似对应这组网格：

```python
GRID_PARAMETERS = {
    "nnnp": 1990,
    "n": 1990,
    "h": 0.015,
    "rnt_scale": 2e-6,
    "hp": 0.0,
    "accy": None,
    "point_n": None,
    "point_h": None,
    "point_rnt_scale": None,
}
```

若严格重现手册的四舍五入 RNT，系数应为 `92*2.17e-8=1.9964e-6`。
该示例交互记录的 ACCY 仍为 `1.5625e-8`。当前补丁会在修改 H 后更新默认
ACCY；要匹配示例记录，须显式填写 `accy=1.5625e-8`。
同时匹配核数据、CSF、ASF 权重、轨道选择和所有运行输入，才能比较其结果。
不能把手册示例与新的自适应 ACCY 设置视为完全相同的数值实验。

## 5. 实际执行的参数试验矩阵

先固定核模型、同位素、CSF、目标态、ASF 权重、变分轨道、优化顺序、阻尼和
物理修正。除被测试项外，不在同一轮改变其他参数。
第一轮可选小模型；通过后，再在最终活性空间复验。

每组完整的 Python 变量设置见第 13 节，可直接替换脚本顶部用户配置区。

### 5.0 B0：未经修改的原版基线

先独立编译并运行原版 B0，**不对它应用补丁**。参考原版的有限核默认值是：

| NNNP | NNN1 | N | H | rnt_scale | HP | ACCY |
| --- | --- | --- | --- | --- | --- | --- |
| 590 | 600 | 590 | 0.05 | `2e-6` | 0 | `H**6=1.5625e-8` |

保存原版源码版本、构建环境和实际输出参数，使用与修改组相同的物理输入。
若原有计算通过交互输入改过网格或 ACCY，另记为历史配置，不冒充默认 B0。
若 B0 未收敛或失败，也应保留记录；它可作为修改前的行为证据，但不是精确解。

这里须区分两种比较：B0 对修改组用于判断相对于原版的整体变化；
下面的 C/D 系列用于在相同补丁处理和固定 ACCY 下分析容量与分辨率。
C0 的网格与 B0 相同，但 C0 已应用补丁并改了 ACCY，因此两者不是同一组。

### 5.1 容量控制与 H 收敛

此表每组都从 B0 的原版源码复制后应用补丁，均采用有限核、
`rnt_scale=2e-6、HP=0、ACCY=1e-10`。表中的 2990 不是原版默认值。

| 编号 | NNNP | N | H | 比较与目的 |
| --- | --- | --- | --- | --- |
| C0 | 590 | 590 | 0.05 | 保留原版容量和网格，应用补丁并固定 ACCY，作为修改组的粗网格对照 |
| C1 / D0 | 2990 | 590 | 0.05 | 对比 C0，仅改变容量；确认容量没有暗中改变实际网格 |
| D1 | 2990 | 1179 | 0.025 | 对比 D0，加密但保持外端 |
| D2 | 2990 | 1965 | 0.015 | 对比 D1，继续加密；外端约增加 1%，应记录这一小差别 |
| D3 | 2990 | 2946 | 0.010 | 对比 D2，验证较细网格平台 |
| D4，按需 | 3990 | 3928 | 0.0075 | 前述变化仍不满足目标时再加密，不预先认为 2990 容量足够 |

B0→C0 同时包含补丁行为调整和 ACCY 调整，不能将差异全部归为网格加密，
也不能将它当作只改变 ACCY 的严格单因素试验。
C0/C1 的输入网格及 ACCY 完全相同，预期结果在数值重复性容差内一致。
D0–D3 属于物理收敛试验，不能要求不同网格的波函数文件逐字节一致。
一旦 C0/C1 的容量对照通过，C0 的结果可作为 D0 的数值对照复用，
但两份程序的容量、资源消耗和源码记录仍分别保存。
至少观察连续两次加密的变化；如果只有最后一对足够小而前一对明显不稳定，
增加一级或采用更小误差预算，不凭一次巧合抵消宣布收敛。

### 5.2 外端范围试验

以 D2 为例，固定 `H=0.015、rnt_scale=2e-6、HP=0、ACCY=1e-10、NNNP=2990`：

| 编号 | N | R(N) / D2 的 R(N) | 原因 |
| --- | --- | --- | --- |
| L0 | 1965 | 1 | 起点 |
| L1 | 2012 | 约 2.024 | 检查延长尾部是否影响结果 |
| L2 | 2058 | 约 4.035 | 再次验证远端截断 |

若采用别的 H，请重新按公式计算 N，不直接套用本表。
要记录最弥散的谱学和关联轨道：MTP、平均半径、尾部概率、节点及关键性质。
R(N) 很大但轨道 MTP 不变时，重点是结果是否稳定和是否增加了无效计算开销。

### 5.3 近核尺度试验

固定 `H=0.015、HP=0、ACCY=1e-10、NNNP=2990`，补偿 N 以维持相近外端：

| 编号 | rnt_scale | N | R(N) / 原 D2 的 R(N) | 原因 |
| --- | --- | --- | --- | --- |
| O0 | `2e-6` | 1965 | 1 | 基准 |
| O1 | `1e-6` | 2012 | 约 1.012 | 减小第一个非零节点及近核尺度 |
| O2 | `5e-7` | 2058 | 约 1.009 | 再减小尺度，验证内层敏感量 |

必须在同一有限核分布下比较。重点测试核附近密度、超精细常数、场移因子和
内层轨道；只看总能量可能遗漏局部积分误差。

### 5.4 ACCY 与初始化路径试验

在通过范围/近核检查的候选网格上，固定全部网格参数，比较
`ACCY=1e-10 → 1e-11 → 1e-12`。达到迭代上限、异常退出或状态漂移的组不能
作为参考结果；不能为了通过而直接放宽节点筛选。

同一网格至少比较两种初始化：

1. 在该网格上重新运行 `rwfnestimate`，使用相同估计方法。
2. 从已收敛参考轨道插值到该网格后重新优化。

若最终结果不同，应先排查局部极值、根选择和未充分收敛。
轨道加载会通过 [LODRWF](../src/lib/lib9290/lodrwf.f90) / INTRPQ 插值到当前
网格，并不自动证明插值后的状态仍是正确解。

### 5.5 点核试验单独进行

固定 `point_rnt_scale=exp(-65/16)、HP=0`，候选 `(point_n, point_h)` 为
`(220,0.0625)`、`(439,0.03125)`、`(877,0.015625)`，均保持 `(N-1)H=13.6875`。
原版点核默认实际点数为 220，容量仍为 590。P0/P1/P2 从干净原版复制，
统一显式设置 `NNNP=990、NNN1=1000、ACCY=1e-10`；990 足以容纳最大的 877 点。
在三个点核修改组中固定容量，避免在加密时同时改变容量。
P0 是打补丁后的粗网格控制组，不是未经修改的原版点核结果。
若需要点核原版基线，在另一计算目录用未经修改的 B0 程序和点核输入运行，
记录默认 `H=0.0625、N=220、ACCY=H**6`，不能复用有限核 B0 的结果。
有限核配置键设为 `None`，不要同时改变有限核参数。
这些键保留的是原版赋值规则；其中有限核 `N=NNNP` 会跟随扩容，P 系列仅测试点核分支。
点核与有限核之间的能量差包含核模型变化，不能当成网格误差。

### 5.6 先执行哪些组

有限核能级研究先跑 B0、C0、D1、D2、D3：保留原版结果，再看固定 ACCY
下粗、中、细网格的变化。C1 可在首次检查扩容补丁时做一次；若省略它，
C0→D1 同时改变容量和网格，不能声称已单独验证容量不影响结果。
候选网格确定后，再做一组外端范围和一组更严格 ACCY 的检查。
仅在这些检查不稳定时继续追加 D4、L2、A2 等组。
超精细、场移等核敏感目标还需要近核尺度检查；点核组只用于点核模型。

原版网格试验统一使用原版程序及同一算法。suite 的集中参数和轨道优化器
另作对照，不能把原版→suite 的全部差异归为网格变化。

## 6. 必须准备哪些测试数据

### 6.1 输入资料清单

每个体系、每个活性空间必须保留：

- `isodata`：Z、同位素 A、有限核分布参数、核质量；做超精细计算还需 I、磁矩、四极矩。
- CSF 文件及生成输入：完整相对论轨道列表、各 J/宇称块、CSF 数量和顺序。
- `rwfnestimate` 的估计方法和输入；使用参考轨道时保留其来源、原网格和校验和。
- RMCDHF 输入：目标 ASF、权重、变分/谱学轨道列表、更新顺序、最大循环数、阻尼及额外覆盖。
- RCI 输入及物理修正选项；性质程序输入；测试跃迁时保留上下态两套输入及双正交设置。

所有网格组使用相同的物理模型。增加 CSF、修改核半径或更换权重必须另开一轮。
旧 RCI `.res` 保存网格信息；网格试验从新的计算目录启动新计算，不能拿旧恢复文件
代替新网格的验证。角积分仅在 CSF/顺序和相关设置相同时才可复用，并记录校验和。

### 6.2 分层数据集和最低要求

| 层次 | 数据 | 需要观察的内容 | 能证明什么 |
| --- | --- | --- | --- |
| 工程冒烟 | 本仓库 `data/test.c、test.w、test.m、isodata` 的小型 Ni 性质样例 | 程序完成、实际网格、g_J、有限输出 | 新配置和性质流水线可运行 |
| 算法基准 | 单电子点核 1s；再选一个有限核内层样例 | 能量、归一化、节点；点核可与匹配常数和模型的 Dirac 解比较 | 排查基础径向求解错误 |
| 中性多电子 | Ni I 与 Cl I，各至少 AS1、AS2 | 基态/低激发态、精细结构、内外层轨道、正交性、收敛记录 | 同时检查开壳层优化与不同轨道尺度 |
| 弥散轨道 | 目标生产体系中最弥散的关联轨道、较高 n 谱学轨道；存在时加入近阈值态 | 平均半径、尾部、MTP、弱跃迁积分 | 外端范围是否足够 |
| 重元素 | 如果研究重中性体系，加入手册 U I 小构型或自己的目标重元素 | 近核密度、内层径向解、能级与目标性质 | 较大 Z 下能否采用该配置 |
| 最终模型 | 实际发表/生产使用的最大活性空间与全部目标态 | 最后两个通过的网格、ACCY 控制、目标性质 | 能否正式选定生产参数 |

Ni、Cl 并不能替代 U 或其他重体系的验证；小模型不能替代最大活性空间验证。
无需第一步就在所有体系跑最大模型：先工程冒烟，再小模型完整扫描，最后生产模型复验。

仓库中的 `data/isodata` 为 Ni 测试夹具，I、磁矩、四极矩等数值含测试设定，
不能直接作为实验核参数使用。Ni/Cl/Fe 大型回归输入通常在仓库外的
`data/rmcdhf_test_data/inputs/` 或 `GRASP_TEST_DATA_ROOT`；脚本不会生成或下载它们。
U I 数据也须从对应手册示例或自己的模型准备，不能把本仓库的小样例当成 U 数据。

### 6.3 每次计算必须保存的输出

| 数据类别 | 保存文件/数值 | 检查目的 |
| --- | --- | --- |
| 实际设置 | 各阶段 stdout、`.sum/.log` 中的 RNT/H/HP/N/ACCY；由其计算 R(N) | 确认交互输入、初始化和后处理没有使用另一组默认值 |
| 计算完成状态 | 退出码、SCF 循环数、各变分轨道 SCNSTY、最后若干轮能量/轨道变化 | 区分严格收敛、迭代上限和非有限数 |
| 态标识 | J/宇称、主组态/LSJ、主 CI 分量、对应根映射 | 不把换根误认为网格变化 |
| 能量 | 总能量 Eh、每个目标态的激发能、精细结构分裂 cm⁻¹ | 总能量和能量差分别收敛 |
| 轨道 | `rwfn.out` 或 `.w`、各轨道能量、P/Q、MTP、节点、归一化、平均半径 | 排查内层欠分辨、尾部截断和错误节点 |
| 轨道关系 | 同一计算的正交残差；不同网格插值后的相位无关重叠 | 优化与插值是否稳定 |
| 超精细 | A/B 及使用的 I、磁矩、四极矩 | 局部性质收敛；避免核因子混淆 |
| 同位素移位 | NMS/SMS 常数、场移因子、最终目标同位素对的位移 | 单项与抵消后的结果都稳定 |
| 跃迁 | 线强、振子强度/跃迁率、长度和速度形式，强线与敏感弱线 | 目标矩阵元稳定；规范差不能单独证明网格收敛 |
| 资源 | wall time、峰值内存、MPI 进程数、线程数、编译器/BLAS | 评估精度与成本，不混用不同资源条件 |

只测试实际需要的性质，但能量与轨道检查每个体系都必须做。
若超精细或同位素位移是最终目标，就必须测试对应量；不能仅凭能量稳定通过。

径向诊断可计算：

```text
轨道归一化：∫(P²+Q²)dr
平均半径：∫r(P²+Q²)dr / ∫(P²+Q²)dr
尾部概率：∫[r_cut,Rmax](P²+Q²)dr，r_cut 在各组使用相同物理半径
```

尾部用两个固定物理 r_cut 检查；不要用“最后 10 个点”比较不同 H。
对于核附近量，也在相同物理核半径附近比较密度和积分。
不同网格的 P/Q 必须先用受控插值映射到共同区间并检验插值误差，再求重叠。
直接按数组下标逐点相减没有物理意义。

## 7. 验收标准：先给出误差预算，再看结果

下表是**初始工程目标示例**，不是 GRASP 官方精度承诺，也不是适用于所有研究的
硬阈值。运行前按项目目标确定自己的绝对/相对容差，并保存在结果配置中。
例如总的激发能误差预算为 1 cm⁻¹，可以暂给网格误差分配 0.1 cm⁻¹；
求解阈值和并行重复性引起的变化再控制在网格预算的 10% 以内。

| 项目 | 初始目标示例 | 解释 |
| --- | --- | --- |
| 实际参数与完成状态 | 参数匹配；无 NaN/Inf；变分轨道通过程序实际收敛判据 | 必须先满足，否则不参与收敛比较 |
| 总能量变化 | `<=1e-6 Eh`，作诊断 | 重体系总能量很大，不能只看其相对误差；该值不替代能量差要求 |
| 目标激发能 | `<=0.1 cm⁻¹` | 随项目目标收紧 |
| 精细结构间隔 | `<=0.01 cm⁻¹` | 先确认相同物理态和正确顺序 |
| g_J | 绝对变化 `<=1e-6` | 零附近也可使用 |
| HFS、场移和质量移位参数 | 相对变化 `<=1e-3`，另设相应单位的绝对容差 | 接近零或存在抵消时，相对值不可靠 |
| 选定强线线强/跃迁率 | 相对变化 `<=1e-3` | 弱线另列绝对容差；不能只检查一个强线 |
| 内部高阶积分归一化/正交残差 | 暂以 `1e-8` 为诊断目标 | 与积分精度、轨道类型相匹配；不能用低阶外部积分直接验收此目标 |
| 节点与态标识 | 保持预期谱学节点、相同 J/宇称与物理态 | 关联轨道要结合其优化约束判断，不机械套用谱学标准 |

能量单位换算可沿用仓库比较工具的 `1 Eh ≈ 219474.6313705 cm⁻¹`。
总能量门槛和激发能门槛需分别检查，不能用一个代替另一个。

使用容差判断时，明确：

```text
|X_new-X_ref| <= max(tol_abs, tol_rel*|X_ref|)
```

对接近零的量必须填写非零、具有物理单位的 `tol_abs`。
最终通过条件是：正确收敛和态匹配；H、范围、近核尺度三项都通过；
ACCY 与 MPI 重复性变化足够小；最后网格组合在生产模型复验通过。
网格偏差不保证单调下降，不采用“能量更低就一定更准”的判定。
与实验数据一致只能作模型核对，不能单独证明径向离散已收敛。

## 8. 现有测试工具怎么使用，哪些不能直接用

### 8.1 脚本与构建检查

本工作区使用共享 Python 环境：

```sh
cd grasp-atomic-suite
../graspkit-tools/.venv/bin/python test/test_patch_grasp_grid.py
```

这些离线测试检查配置入口、重复定义、范围冲突、备份、幂等和恢复，不进行
原子物理计算。CMake/CTest 通过也不等于目标体系已达到网格收敛。

### 8.2 能量比较

对于相同物理态、相同 CSF 和可由工具解析的 RMCDHF 输出，可在服务器运行：

```fish
python3 /path/to/grasp-atomic-suite/test/rmcdhf_orbopt/compare_sum.py \
    /path/to/results/D3/rmcdhf.sum \
    /path/to/results/D2/rmcdhf.sum \
    --allow-radial-grid-difference \
    --tolerance 1e-6
```

该命令只是本文示例总能量阈值的检查；应另算激发能、精细结构和最终性质。
只有确认根匹配后才可比较相同序号；工具不会自动识别全部态交换。
`--allow-energy-differences` 仅适合探索时输出差值，不能作为最终通过条件。

`compare_fine_structure.py` 可比较特定 J/宇称块的 level-1 间隔：

```fish
python3 /path/to/grasp-atomic-suite/test/rmcdhf_orbopt/compare_fine_structure.py \
    D2=/path/to/results/D2/rmcdhf.sum \
    D3=/path/to/results/D3/rmcdhf.sum \
    --j-values 2,3,4 --parity +
```

这里只是对应 Ni 目标序列的用法示例。其他体系按真实 J/宇称与根修改；
这个工具不自动处理每块内多个目标根，也不自动执行物理验收。

### 8.3 波函数与集群回归的限制

`compare_rwfn.py` 的重叠计算要求相同节点前缀，跨 H/RNT 网格时会拒绝。
它采用梯形积分，适合独立趋势检查，不是 GRASP QUAD 的高精度替代。
跨网格收敛需受控插值和积分工具；当前仓库没有覆盖本文全部指标的一键验收程序。
不能把该脚本的低阶归一化偏差直接与本节的 `1e-8` 内部积分目标比较。

`run_data_case.sh` 和 SLURM 模板用于现有 orbopt 回归，部分脚本会加载指定的
GRASP module，并从 PATH 找 `rwfnestimate/rangular_mpi`。
它们不能无修改地当成任意原版 GRASP 副本的一键网格测试器。
使用时必须核对 module、所有程序绝对路径、外部输入和初始轨道；只设置
`GRASP_BINDIR` 不保证准备程序也来自同一份原版。
原版 `rmcdhf_mpi` 与 suite `rmcdhf_orbopt_mpi` 也不能直接混用名称。

## 9. fish 下每份源码的操作流程

保留 B0，并从同一份未经修改的原版源码分别复制其他组，例如：

```text
grasp-B0/  # 原版 590/600，不应用补丁
grasp-C0/
grasp-D0/
grasp-D1/
grasp-D2/
grasp-D3/
results/<体系>/<活性空间>/<试验编号>/
```

复制时不带入旧 `build/`、编译产物或源码内 `.mod/.o`；CMake 缓存含绝对路径。
不要从当前已调整到 2990 的源码、suite 或上一个已打补丁的测试组继续复制。
若手头只有调整后的副本，先从确认的原版提交或原始发布包取得干净源码。
修改组的每份根目录保存一个 `patch_grasp_grid.py` 副本，便于参数与程序对应。
B0 跳过补丁步骤，直接按下面相同的 CMake 流程编译。
下面是 D2 的完整文件配置：

```python
# D2：第 9 节操作示例
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

在 fish 中执行，不传数值参数：

```fish
cd /path/to/grasp-D2
python3 patch_grasp_grid.py
```

查看预览后，将文件中的 `RUN_MODE` 改成 `"apply"`，执行同一命令；
再改成 `"check"`，执行同一命令，确认 `Files requiring changes: 0`。
备份自动保存在该份源码的 `grid-backups/<时间戳>/`；同时保留此次配置和预览输出。

然后加载服务器 Fortran/MPI/BLAS/LAPACK 环境，按
[原版 README 的 CMake 流程](https://github.com/compas/grasp/blob/master/README.md#cmake-based-build)：

```fish
cd /path/to/grasp-D2
./configure.sh
cd build
make -j4 install
ctest --output-on-failure
```

`./configure.sh` 只用于尚无 `build/` 的首次配置；同一份以后改参数时：

```fish
cd /path/to/grasp-D2/build
make clean
make -j4 install
ctest --output-on-failure
```

MPI 没有被 CMake 找到时不能宣称 MPI 版本已经更新。
编译完成后，在独立计算目录通过 `/path/to/grasp-D2/bin/<程序名>` 调用。
至少核对 `rwfnestimate → rmcdhf[_mpi] → rci[_mpi] → 目标性质` 所有阶段
的程序来源及实际设置，不能只核对最后一个程序。

suite 模式将 `SOURCE_LAYOUT` 改为 `"atomic-suite"`，指向 suite 根目录，
脚本只调整公共文件；轨道命令是 `rmcdhf_orbopt*`，构建需启用实际使用的
串行/MPI 目标。网格数值和物理验证规则与原版相同。

## 10. 结果记录模板和最终选型

每个试验保存脚本配置、源码提交/差异、二进制 SHA256、编译配置、module 清单、
输入文件校验和和完整 stdout。建议有两个 CSV，而不是把不同物理态挤在一行：

```text
# runs.csv：每次运行一行
run_id,dataset,isotope,active_space,nnnp,n,h,hp,rnt_scale,actual_rnt,actual_accy,rmax_a0,nprocs,threads,exit_code,strict_scf_converged,iterations,wall_seconds,peak_memory,source_revision,binary_sha256

# observables.csv：每个目标量一行
run_id,quantity,state_or_transition_id,value,unit,reference_run,delta_abs,delta_rel,tol_abs,tol_rel,pass
```

物理态标识包括 J/宇称、项或主组态和对应根，不能只写“第 3 个能级”。
`pass` 应依据事先约定的容差生成；失败项保留，不覆盖旧结果。

最后选**满足全部目标量、验证了尾部与初始化、成本可接受**的最小网格。
报告同时列出 B0 原版结果、相对于 B0 的变化，以及修改组之间的收敛变化。
与 B0 一致本身不能证明准确；若 B0 网格不足，改进后的结果可以偏离 B0。
如 D2 与 D3、范围试验、近核试验、ACCY 试验及生产模型均通过，可选 D2；
否则按最敏感性质继续加密或扩容。若单项优选后形成新的 H/N/RNT 组合，
再做一次组合复验，因为这些误差不保证相互独立。

串行、MPI 1/2/4 进程以及实际使用的内存版本，至少在小模型上做同网格一致性
检查；生产进程规模另做一次复验。允许浮点求和和迭代顺序导致的小变化，
但必须低于重复性预算并保持相同物理态，不能用“MPI 不必一致”跳过异常。

## 11. 当前已有证据与尚需完成的工作

- 已有：20 个脚本离线用例通过；两份真实原版源码副本的无参数配置、独立应用、
  检查、重复应用和逐字节恢复通过；这些验证未修改参考原版源码。
- 已有：suite 10 个应用的构建、11 项 MPI 构建 CTest、10 项串行 CTest，以及
  小型 Ni g_J 样例的历史一致性验证。历史记录见
  [公共参数说明](common_parameters.md)与[迁移说明](rmcdhf_migration.md)。
- **尚未完成**：本文 D/L/O/ACCY 矩阵对应的 Ni/Cl/U/生产体系物理计算、
  全部目标性质的误差预算验收及目标服务器上的 MPI 规模复验。

只有完成后者，才可以在结果报告中写“该体系的网格已达到指定精度”。

## 12. 依据与相关文档

- 核查基准：本工作区参考原版 `parameter_def_M.f90`、`getscd.f90`、
  `radgrd.f90`、`scf.f90`、`count.f90`、`lodrwf.f90`；本文公式和试验组合据此推导。
- 原版手册：[Jönsson 等，GRASP Manual for Users，Atoms 2023, 11, 68](https://www.diva-portal.org/smash/get/diva2:1771324/FULLTEXT01.pdf)，
  §1.4、§13.4、§13.5。手册给出扩容和具体实例，没有为本文验收阈值作保证。
- 构建：[compas/grasp README](https://github.com/compas/grasp/blob/master/README.md#cmake-based-build)。
- 接口与备份：[脚本操作说明](../scripts/README.md)。
- suite 配置：[公共参数说明](common_parameters.md)。
- 回归与历史限制：[轨道测试说明](../test/rmcdhf_orbopt/README.md)、
  [历史测试结果](../test/rmcdhf_orbopt/RESULTS.md)。


## 13. 可直接复制的 Python 配置

### 13.1 复制规则与索引

以下每个 Python 块都是脚本顶部用户配置的**完整替换内容**，不是新的配置文件。
第 13.3 节 B0 的 `ORIGINAL_BASELINE` 是原版参数记录，单独标明，不属于补丁配置。
只选一组，替换 `GRASP_SOURCE` 到 `RESTORE_DIRECTORY` 的配置区，不修改下面的实现。
不要把多组代码依次粘进同一个脚本，否则后面的赋值会覆盖前面的配置。

每组均假定脚本已经复制到对应 GRASP 根目录，因此 `GRASP_SOURCE="."`。
如果脚本仍放在 suite 的 `scripts/` 中，就把该值替换为目标源码的绝对路径。
原版各组均从同一份未经修改、容量为 590 的 GRASP 复制，使用独立计算结果目录。
只有第 13.9 节的 suite 示例从 suite 源码复制，不能替代原版 B0。

这些块默认 `RUN_MODE="preview"`。运行后改成 `"apply"`，再次运行写入；
然后改成 `"check"` 运行检查。不要重新粘贴块导致模式又回到 preview 而误以为已写入。
各组 `RESTORE_DIRECTORY=None`，避免继承上一组的恢复设置。
完成修改后，按第 9 节手动编译；复制参数本身不会更新二进制。

| 测试 | 配置块所在小节 | 对照关系 |
| --- | --- | --- |
| 原版基线 | 13.3 B0 | 未经修改的原版，先编译运行；不应用补丁 |
| 手册对照 | 13.2 M-U | 单独匹配手册输入，不混入 D 系列 |
| 容量 | 13.3 C0、C1 | C0 ↔ C1 |
| 分辨率 | 13.3 D0–D4 | D0 → D1 → D2 → D3，按需 D4 |
| 外端范围 | 13.4 L0–L2 | L0 → L1 → L2 |
| 近核尺度 | 13.5 O0–O2 | O0 → O1 → O2 |
| 数值阈值 | 13.6 A0–A2 | A0 → A1 → A2 |
| 点核分辨率 | 13.7 P0–P2 | P0 → P1 → P2 |
| 初始化 | 13.8 I-EST、I-REF | 同网格，不同初始化来源 |
| 串行/MPI | 13.8 R 系列 | 同网格、同模型，不同程序或进程数 |
| suite 使用 | 13.9 | 相同数值，切换源码布局 |

下面列出**每组最终填写的值**，以及实际比较时改变的参数。所有原版修改组
都独立从原版 590 点源码复制，而不是在前一组源码上继续修改。
有限核 D0–D3/L/O/A/I/R 系列统一容量 2990，按需 D4 单独扩容到 3990；
P 系列统一容量 990；
二者都属于显式扩容后的测试设置。脚本自动把 `NNN1` 设为容量加 10。

| 组别 | nnnp → NNN1 | 实际点数 | H | RNT 系数 | ACCY | 比较对象及改变项 |
| --- | --- | --- | --- | --- | --- | --- |
| B0 | 590 → 600 | 590 | 0.05 | `2e-6` | `1.5625e-8` | 未经修改的原版；仅记录，不应用补丁 |
| C0 | 590 → 600 | 590 | 0.05 | `2e-6` | `1e-10` | 相对 B0：补丁处理及 ACCY；保留原版容量与网格 |
| C1 / D0 | 2990 → 3000 | 590 | 0.05 | `2e-6` | `1e-10` | 相对 C0：只扩容，实际网格不变 |
| D1 | 2990 → 3000 | 1179 | 0.025 | `2e-6` | `1e-10` | 相对 D0：N/H 配套加密 |
| D2 | 2990 → 3000 | 1965 | 0.015 | `2e-6` | `1e-10` | 相对 D1：N/H 配套加密 |
| D3 | 2990 → 3000 | 2946 | 0.010 | `2e-6` | `1e-10` | 相对 D2：N/H 配套加密 |
| D4，按需 | 3990 → 4000 | 3928 | 0.0075 | `2e-6` | `1e-10` | 相对 D3：再扩容并加密，包含容量变化 |
| L0 / O0 / A0 / I / R | 2990 → 3000 | 1965 | 0.015 | `2e-6` | `1e-10` | 复用 D2 网格；I/R 另改初始化或运行方式 |
| L1 | 2990 → 3000 | 2012 | 0.015 | `2e-6` | `1e-10` | 相对 L0：只增大 N，延长外端 |
| L2 | 2990 → 3000 | 2058 | 0.015 | `2e-6` | `1e-10` | 相对 L1：继续延长外端 |
| O1 | 2990 → 3000 | 2012 | 0.015 | `1e-6` | `1e-10` | 相对 O0：减半 RNT 系数并补偿 N |
| O2 | 2990 → 3000 | 2058 | 0.015 | `5e-7` | `1e-10` | 相对 O1：再减半 RNT 系数并补偿 N |
| A1 | 2990 → 3000 | 1965 | 0.015 | `2e-6` | `1e-11` | 相对 A0：只收紧 ACCY |
| A2 | 2990 → 3000 | 1965 | 0.015 | `2e-6` | `1e-12` | 相对 A1：只收紧 ACCY |
| P0 | 990 → 1000 | 点核 220 | 0.0625 | `exp(-65/16)` | `1e-10` | 原版点核网格，另有扩容、补丁及 ACCY 调整 |
| P1 | 990 → 1000 | 点核 439 | 0.03125 | `exp(-65/16)` | `1e-10` | 相对 P0：只配套改变 point_n/point_h |
| P2 | 990 → 1000 | 点核 877 | 0.015625 | `exp(-65/16)` | `1e-10` | 相对 P1：只配套改变 point_n/point_h |
| M-U | 1990 → 2000 | 1990 | 0.015 | `1.9964e-6` | `1.5625e-8` | 单独对照手册 U I 记录 |

各组 HP 均为 0。有限核各块的 `point_* = None` 保留原版点核规则，
而点核各块的有限核键为 `None`。每个配置块中的注释进一步标明相对于原版的修改；
未改变数值也显式写出，方便核对该副本最终应使用的设置。

C1 与 D0 相同；D2、L0、O0、A0 和下面 I/R 的示例数值相同。
为方便分别复制，L0/O0/A0 仍各自列出完整代码。对于同一体系和相同运行条件，
可以复用同一已验证基准结果，但必须记录对应关系，不把它计为多个独立结果。

### 13.2 手册 U I 参数对照

此组用于对照第 4 节的手册参数记录，不加入 D 系列的控制变量比较。固定 U I、Z=92、A=238 和同一物理模型；这里显式匹配手册记录的 RNT 与 ACCY。

#### M-U：手册参数对照

检查实际 RNT=2.17e-8、N=1990、H=0.015，以及手册对应的态和输入条件。

```python
# M-U：手册参数对照
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 1990,  # 原版 590 → 1990，扩容；NNN1 自动为 2000
    "n": 1990,  # 原版 590 → 1990，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 1.9964e-6,  # 原值 2e-6 → 1.9964e-06；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1.5625e-8,  # 原版默认 H**6 → 固定 1.5625e-08
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.3 容量与步长测试

先运行 B0。C0 保留 590 点网格但应用补丁并固定 ACCY，不等同于 B0。
C0 与 C1 只改变容量；C1 同时作为 D0。D1→D2→D3 改变分辨率；D4 仅在需要时执行。

#### B0：未经修改的原版，590 点容量

**不用 `patch_grasp_grid.py` 写入这一组。** 直接编译干净原版，用其默认值运行。
以下仅供记录和核对输出，不复制到脚本配置区：

```python
ORIGINAL_BASELINE = {
    "nnnp": 590,
    "nnn1": 600,
    "n": 590,
    "h": 0.05,
    "rnt_scale": 2e-6,
    "hp": 0.0,
    "accy": 1.5625e-8,  # 原版默认 H**6
}
```

将补丁参数填回这些数值再执行 apply，仍不等于保留原版源码：
脚本还会恢复 `rwfnrelabel` 的初始化、调整 RMCDHF 的 ACCY 更新位置等。
`None` 只表示保留对应参数，不能用来把一次 apply 当作未经修改的原版。

#### C0：590 点容量控制组

这是从 B0 复制并打补丁的控制组。先记录相对于 B0 的变化，再对照 C1，
检查相同实际网格和 ACCY 下的能量、轨道与性质是否一致。

```python
# C0：590 点容量控制组
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 590,  # 原版 590，保留容量；NNN1 自动为 600
    "n": 590,  # 原版 590，保留实际点数
    "h": 0.05,  # 原值 0.05，保留步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### C1 / D0：扩容但保持实际网格

与 C0 比较容量影响；随后作为 D1 的粗网格对照。

```python
# C1 / D0：扩容但保持实际网格
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 590,  # 原版 590，保留实际点数
    "h": 0.05,  # 原值 0.05，保留步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### D1：H=0.025

对照 D0，保存目标能级、精细结构、轨道与所需性质。

```python
# D1：H=0.025
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1179,  # 原版 590 → 1179，实际点数
    "h": 0.025,  # 原值 0.05 → 0.025，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### D2：H=0.015

对照 D1，并作为下面范围、近核、ACCY 示例的基准。

```python
# D2：H=0.015
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### D3：H=0.010

对照 D2，检查最后加密的变化是否低于误差预算。

```python
# D3：H=0.010
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 2946,  # 原版 590 → 2946，实际点数
    "h": 0.010,  # 原值 0.05 → 0.01，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### D4：按需进一步加密

对照 D3；必须使用匹配的新容量重新编译全部程序。

```python
# D4：按需进一步加密
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 3990,  # 原版 590 → 3990，扩容；NNN1 自动为 4000
    "n": 3928,  # 原版 590 → 3928，实际点数
    "h": 0.0075,  # 原值 0.05 → 0.0075，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.4 外端范围测试

比较 L0→L1→L2。H、RNT、HP 和 ACCY 保持相同，只延长外端。重点看最弥散轨道的 MTP、相同物理半径处的尾部概率，以及目标物理量。

#### L0：范围基准（数值等同 D2）

作为 L1/L2 的对照，不改变 D2 的初始轨道与模型。

```python
# L0：范围基准（数值等同 D2）
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### L1：约两倍外端

与 L0 比较，记录轨道尾部及性质变化。

```python
# L1：约两倍外端
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 2012,  # 原版 590 → 2012，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### L2：约四倍外端

与 L1 和 L0 比较，验证截断变化是否进一步减小。

```python
# L2：约四倍外端
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 2058,  # 原版 590 → 2058，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.5 近核尺度测试

比较 O0→O1→O2。减小 rnt_scale，同时补偿 N。各组仍使用同一有限核分布；重点看核附近密度、内层轨道、超精细常数和场移因子。

#### O0：近核基准（数值等同 D2）

作为 O1/O2 的对照。

```python
# O0：近核基准（数值等同 D2）
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### O1：RNT 系数减半

与 O0 比较；外端约增加 1.2%，记录其余量是否稳定。

```python
# O1：RNT 系数减半
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 2012,  # 原版 590 → 2012，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 1e-6,  # 原值 2e-6 → 1e-06；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### O2：RNT 系数减为四分之一

与 O1 比较；外端约增加 0.9%，核敏感量必须一起检查。

```python
# O2：RNT 系数减为四分之一
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 2058,  # 原版 590 → 2058，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 5e-7,  # 原值 2e-6 → 5e-07；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.6 ACCY 测试

下面固定 D2 网格作示例，比较 A0→A1→A2。若最终候选是其他网格，要在三组中同时替换 nnnp/n/h/rnt_scale，保证它们除 accy 外完全相同。

#### A0：ACCY=1e-10（数值等同 D2）

记录严格收敛状态、迭代次数、目标物理量。

```python
# A0：ACCY=1e-10（数值等同 D2）
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### A1：ACCY=1e-11

对照 A0，检查求解变化是否小于分配给它的误差预算。

```python
# A1：ACCY=1e-11
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-11,  # 原版默认 H**6 → 固定 1e-11
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### A2：ACCY=1e-12

对照 A1；若不收敛，不以此组作参考，也不直接放宽节点判据。

```python
# A2：ACCY=1e-12
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-12,  # 原版默认 H**6 → 固定 1e-12
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.7 点核测试

仅在 isodata 选择点核模型时使用。从原版 590 点容量源码复制，三组均扩容到
990（NNN1=1000），足以容纳最大的 877 点，保持比较过程中容量不变。
有限核键 n/h/rnt_scale 均为 None；明确固定点核 RNT 系数为 exp(-65/16)
的双精度近似值。比较 P0→P1→P2，保持相同核模型和目标态，
不能与有限核 D 系列直接作为网格误差比较。

#### P0：点核修改组的粗网格控制

点核 N/H 保留原版值，但显式扩容并固定 ACCY。检查能量、节点和轨道积分。

```python
# P0：点核修改组的粗网格控制
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 990,  # 原版 590 → 990，扩容；NNN1 自动为 1000
    "n": None,  # 保留有限核 N=NNNP 规则；本组仅使用点核
    "h": None,  # 保留有限核 H；本组仅使用点核
    "rnt_scale": None,  # 保留有限核 RNT 系数；本组仅使用点核
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版点核默认 H**6 → 固定 1e-10
    "point_n": 220,  # 原版实际点数 220，保留
    "point_h": 0.0625,  # 原版步长 0.0625，保留
    "point_rnt_scale": 0.017205950425851383,  # 原版 exp(-65/16) 的数值表示；RNT=此值/Z
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### P1：点核步长减半

对照 P0；实际点核范围相同。

```python
# P1：点核步长减半
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 990,  # 原版 590 → 990，扩容；NNN1 自动为 1000
    "n": None,  # 保留有限核 N=NNNP 规则；本组仅使用点核
    "h": None,  # 保留有限核 H；本组仅使用点核
    "rnt_scale": None,  # 保留有限核 RNT 系数；本组仅使用点核
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版点核默认 H**6 → 固定 1e-10
    "point_n": 439,  # 原版实际点数 220 → 439
    "point_h": 0.03125,  # 原版步长 0.0625 → 0.03125
    "point_rnt_scale": 0.017205950425851383,  # 原版 exp(-65/16) 的数值表示；RNT=此值/Z
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### P2：点核再次加密

对照 P1；另行验证 ACCY，不把核模型差当作网格误差。

```python
# P2：点核再次加密
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 990,  # 原版 590 → 990，扩容；NNN1 自动为 1000
    "n": None,  # 保留有限核 N=NNNP 规则；本组仅使用点核
    "h": None,  # 保留有限核 H；本组仅使用点核
    "rnt_scale": None,  # 保留有限核 RNT 系数；本组仅使用点核
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版点核默认 H**6 → 固定 1e-10
    "point_n": 877,  # 原版实际点数 220 → 877
    "point_h": 0.015625,  # 原版步长 0.0625 → 0.015625
    "point_rnt_scale": 0.017205950425851383,  # 原版 exp(-65/16) 的数值表示；RNT=此值/Z
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.8 初始化与并行一致性测试

这两类测试不改变网格，所以同一套配置用于所有对照计算；需要改变的是运行输入或程序/进程数。

#### I-EST / I-REF：初始化路径的共用配置

I-EST 在本组网格重新生成估计；I-REF 读入参考轨道并插值后重新优化。两组仅改变初始化来源，保持相同物理模型、变分轨道和求解设置，比较最终物理态和物理量。

```python
# I-EST / I-REF：初始化路径的共用配置
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

#### R-SERIAL / R-MPI1 / R-MPI2 / R-MPI4：并行测试的共用配置

R-SERIAL 使用串行程序；其余使用相同包的 MPI 程序，分别设 1/2/4 进程。使用相同输入与网格，不在 Python 参数中添加 MPI 进程字段。

```python
# R-SERIAL / R-MPI1 / R-MPI2 / R-MPI4：并行测试的共用配置
GRASP_SOURCE = "."
SOURCE_LAYOUT = "grasp2018"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # 原版 590 → 2990，扩容；NNN1 自动为 3000
    "n": 1965,  # 原版 590 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # 原版默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

### 13.9 在 grasp-atomic-suite 中使用同一组数值

上面的修改组可将 `SOURCE_LAYOUT="grasp2018"` 改为 `"atomic-suite"`，并让
`GRASP_SOURCE` 指向那份 suite 根目录；B0 的原版记录不适用此替换。
每份 suite 副本从同一 suite 版本复制，不能把原版 590 点基底的说明照搬为
suite 的默认值。下面是 suite 的 D2 完整示例：

```python
# SUITE-D2：suite 中采用 D2 网格
GRASP_SOURCE = "."
SOURCE_LAYOUT = "atomic-suite"
RUN_MODE = "preview"
PRINT_DIFF = True
GRID_PARAMETERS = {
    "nnnp": 2990,  # suite 当前默认 2990，保留容量；NNN1 自动为 3000
    "n": 1965,  # suite 当前默认 2990 → 1965，实际点数
    "h": 0.015,  # 原值 0.05 → 0.015，减小步长
    "rnt_scale": 2e-6,  # 原值 2e-6，保留；RNT=此值/Z
    "hp": 0.0,  # 原值 0，保持指数网格
    "accy": 1e-10,  # suite 当前默认默认 H**6 → 固定 1e-10
    "point_n": None,  # 保留原版 MIN(220,NNNP) 规则
    "point_h": None,  # 保留点核默认 H=0.0625
    "point_rnt_scale": None,  # 保留点核默认 exp(-65/16)
}
BACKUP_DIRECTORY = None
RESTORE_DIRECTORY = None
```

suite 需要使用 `rmcdhf_orbopt*` 等对应程序，并使原版准备/RCI 程序具有匹配
容量和网格。程序选择在作业脚本中设置，不是给 GRID_PARAMETERS 添加字段。
如果采用第 8 节提到的集群模板，还要检查它是否重新加载了另一份 GRASP module。

节点阈值诊断是另一项试验：Python 脚本没有 THRESH/NODE_THRESHOLD 配置键。
若单独改变 suite 的 NODE_THRESHOLD，只能改公共 Fortran 文件并重编译；
在该诊断期间保持所选 GRID_PARAMETERS 完全相同。HP 非零网格也不在本轮
已定义试验矩阵中，不能把 HP 从 0 改成任意正值后直接复用 L/O 表中的范围。

### 13.10 复制后逐项确认

1. 文件只保留一个用户配置块，源码路径和布局正确。
2. 保留并运行未经修改的 B0；原版修改组从同一份 590 点源码复制。
   B0/C/D/L/O/A 系列用同一有限核 isodata；P 系列用点核 isodata。
3. 各组值与第 5 节的矩阵一致，固定 ACCY 的组没有被交互输入覆盖。
4. 实际应用、检查通过并全量编译安装；计算调用该组程序的绝对路径。
5. 保存实际 RNT/H/HP/N/ACCY、目标态、收敛状态与第 6 节的观测数据。
6. 依据第 7 节的误差预算验收。以上配置是测试候选，不是已验证的生产默认值。
