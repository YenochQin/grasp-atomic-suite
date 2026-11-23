# gj / gj_mpi：Landé 因子

`gj` 和 `gj_mpi` 计算选定原子能级的电子 Landé 因子 `g_J`，并给出电子反常磁矩修正。
`gj` 为串行程序，`gj_mpi` 将矩阵元计算分配给多个 MPI 进程。

程序读取 `isodata`、`name_desc.parquet`、`name_header.toml`、`name.wfn` 和混合系数文件。
MCDHF 模式读取 `name.mix` 并写出 `name.gj`，CI 模式读取 `name.cimix` 并写出 `name.cgj`。
`gj_mpi` 支持交互输入，也可从命令行接收状态名及 `--nonci` 或 `--ci`。

[RHFS 程序](../rhfs/README.md) · [Landé 因子计算流程](../../../docs/implemented/design-0027-rhfs-lande-g.md)

## 计算流程

1. 从混合文件取得选定能级的 CSF 系数 $c_i$。在共同径向网格上，对每个轨道对预制磁矩径向积分和径向重叠，并结合一电子角因子生成普通磁矩与反常磁矩两张轨道表。[输入与串行入口](gj.f90)
2. 对每对同 $J$ 的 CSF 计算一电子张量角系数，将它们与轨道表相乘，得到两套 CSF 矩阵元。对角项用 $c_i^2$，下三角非对角项合并 $c_ic_j+c_jc_i$ 后累加到各态。MPI 版本按列分工并分别归约两张结果数组。[串行累加](compute_gj_only.f90) · [MPI 累加](compute_gj_only_mpi.f90)
3. 对 $J>0$ 的能级，分别将普通和修正约化矩阵元乘 $c_0/\sqrt{J(J+1)}$ 和 $0.001160/\sqrt{J(J+1)}$，写出两项及其和。这里的 $c_0$ 为原子单位下的真空光速。[归一化与报告](write_lande_g_factors.f90)

所需的径向超精细积分、单电子角因子和部分输入逻辑与 RHFS 共用，目标的具体源文件由构建清单指定。[CMake 构建清单](CMakeLists.txt) · [RHFS 各文件流程](../rhfs/README.md)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`gj.f90`](gj.f90) | 交互选择状态与 MCDHF／CI 混合文件，加载 CSF、核和径向轨道。计算所选本征态的两个磁矩期望值并输出。 |
| [`gj_mpi.f90`](gj_mpi.f90) | rank 0 解析命令行或交互输入并广播状态，集体载入 CSF。各 rank 读取径向和混合数据，归约两个期望值，由 rank 0 写报告。 |

### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`compute_gj_only.f90`](compute_gj_only.f90) | 预先计算轨道对的超精细径向积分、重叠与磁矩角因子。扫过下三角同 $J$ 的 CSF 对，用单电子张量系数组装磁矩和自旋修正，再按混合系数乘积累加到每个能级的对角元。 |
| [`compute_gj_only_mpi.f90`](compute_gj_only_mpi.f90) | 沿 CSF 列循环分派任务，各 rank 预制轨道表并分别累加态上的磁矩、反常磁矩项。两条数组分别作 MPI 全归约。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`write_lande_g_factors.f90`](write_lande_g_factors.f90) | 由 `2J+1` 还原 $J$，跳过 $J=0$。用 $1/\sqrt{J(J+1)}$ 归一化两类约化矩阵元，分别乘原子单位真空光速和 `0.001160`，写出 $g_J$、修正和总和。 |
