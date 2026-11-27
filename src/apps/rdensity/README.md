# rdensity / rdensity_mpi：径向电子密度与自然轨道

`rdensity` 和 `rdensity_mpi` 计算各选定能级的径向电子密度，并对这些能级的单电子密度矩阵
按统计权重平均，得到自然轨道。`rdensity_mpi` 的密度矩阵计算由多个 MPI 进程分担。

程序读取 `isodata`、`name_desc.parquet`、`name_header.toml`、`name.wfn` 和
`name.mix` 或 `name.cimix`，将自然轨道写入 `name.mcdhf.natural-orbitals`（MCDHF）或 `name.ci.natural-orbitals`（CI），将径向密度写入
`name.density`（MCDHF）或 `name.cidensity`（CI）。MPI 模式由 rank 0 写出这些文件。

[输入选项与自然轨道排序](../../../docs/implemented/design-0026-rdensity-mpi.md)

## 计算流程

1. 读入径向波函数、原生 CSF 基组及每个原子能级的混合系数。按同 $J$、宇称以及一电子算符允许的占据差筛选 CSF 对，收缩为各态的同 $\kappa$ 一体密度矩阵 $D^{(s)}$。MPI 路径让各 rank 计算不同列，并将各级密度块求和到 rank 0。[密度组装](density_calculation.f90) · [MPI 分工](density_parallel_mpi.f90)
2. 先用原始轨道收缩 $D(r)=\sum_{ab}D^{(s)}_{ab}(P_aP_b+Q_aQ_b)$，写出正值网格点上的 $D(r)$ 与 $\rho(r)=D(r)/r^2$。再按 $2J_s+1$ 对各级密度矩阵做归一化加权平均。[密度与统计平均](density_calculation.f90)
3. 每个 $\kappa$ 块调用 LAPACK `DSYEV`，得到自然占据数及轨道旋转矩阵。按指定次序重排本征向量，以同一旋转处理大、小分量及原点系数，并写出 `RAD--WFN` 径向流。[对角化与旋转](natural_orbitals.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 初始化诊断位，读取允许的非默认诊断选择，按需连接 `.dbg` 单元，供输入、积分与结果检查使用。 |
| [`rdensity.f90`](rdensity.f90) | 解析状态、混合文件来源及自然轨道排序方式，检查输入并装入基组和径向轨道。初始化后调用密度组装，MPI 与串行都用同一入口。 |
| [`read_mixing_blocks.f90`](read_mixing_blocks.f90) | 读取状态对应的 `CSF--MIX`，逐块提取根编号、$J$、宇称、绝对能量和本征向量，将其按 CSF 偏移铺到共享数组，供各算符按态收缩。 |
| [`read_radial_density_input.f90`](read_radial_density_input.f90) | 读取核数据、光速和径向网格设置，初始化径向差分和网格，载入 `name.wfn` 中用于构造密度的径向函数。 |

### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`density_calculation.f90`](density_calculation.f90) | 按 $J$、宇称和一体占据选择筛选 CSF 对，将角系数与混合乘积累加为同 $\kappa$ 的密度块。先用原轨道收缩得 $D(r)=r^2\rho(r)$，对正值网格点写 $D(r)$ 和 $\rho(r)$。再按 $2J+1$ 平均密度块并旋转到自然轨道，保存为自然轨道径向流。 |
| [`natural_orbitals.f90`](natural_orbitals.f90) | 用 LAPACK `DSYEV` 对每个密度块求本征占据数与旋转矩阵。按选项以占据降序或原轨道主导成分匹配排序。用同一旋转作用于 $P,Q$ 和原点系数。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`density_parallel_mpi.f90`](density_parallel_mpi.f90) | 初始化 MPI 并广播状态、混合来源、排序和错误状态，利用集体 reader 载入 CSF，将每级密度块求和到 rank 0。 |
| [`density_parallel_serial.f90`](density_parallel_serial.f90) | 安装串行 CSF 基组、提供固定 rank 上下文。广播和密度求和没有跨进程数据传输，错误直接返回终止状态。 |
