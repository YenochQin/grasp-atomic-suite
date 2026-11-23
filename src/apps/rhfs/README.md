# rhfs / rhfs_mpi：相对论超精细结构

`rhfs` 和 `rhfs_mpi` 从 MCDHF 或 CI 波函数计算磁偶极超精细常数 `A`、电四极超精细常数
`B`、非对角超精细矩阵元和电子 Landé 因子 `g_J`。`rhfs_mpi` 将矩阵元计算分配给多个 MPI 进程。

程序读取 `isodata`、`name_desc.parquet`、`name_header.toml`、`name.wfn` 和
`name.mix` 或 `name.cimix`。MCDHF 模式写出 `name.h` 和 `name.hoffd`，
CI 模式写出 `name.ch` 和 `name.choffd`。

[独立 Landé 因子程序](../gj/README.md) · [超精细常数与 Landé 因子计算流程](../../../docs/implemented/design-0027-rhfs-lande-g.md)

## 计算流程

1. 装入同位素的核矩与核自旋数据、原生基组、径向轨道和能级混合系数。以轨道对为单位预制一电子角因子和所需径向积分。[输入与配置](read_hyperfine_options.f90) · [角因子](hyperfine_one_electron_element.f90) · [径向积分](hyperfine_radial_integral.f90)
2. 分别对磁偶极、电四极和 Landé 项计算 CSF 约化矩阵元，按初、末能级的混合向量收缩为 ASF 矩阵。乘入核磁偶极矩、核电四极矩及相应的角动量归一化因子，输出对角常数和允许的非对角矩阵元。[串行超精细计算](calculate_hyperfine_matrix_elements.f90) · [超精细与 Landé 归一化说明](../../../docs/implemented/design-0027-rhfs-lande-g.md)
3. MPI 版本由各 rank 负责不同的 CSF 列，归约超精细、Landé 和自旋修正数组后由 rank 0 完成级结果输出。[MPI 超精细计算](calculate_hyperfine_matrix_elements_mpi.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 初始化诊断位，读取允许的非默认诊断选择，按需连接 `.dbg` 单元，供输入、积分与结果检查使用。 |
| [`read_hyperfine_options.f90`](read_hyperfine_options.f90) | 读核自旋和核矩等同位素数据，以及可选光速、网格和一阶 CSF 截止。初始化网格、差分和核势，再把径向轨道装入内存。 |
| [`read_mixing_blocks.f90`](read_mixing_blocks.f90) | 读取状态对应的 `G92MIX`，逐块提取根编号、$J$、宇称、偏移能量和本征向量，将其按 CSF 偏移铺到共享数组，供各算符按态收缩。 |
| [`rhfs_main.f90`](rhfs_main.f90) | 读取状态与混合来源，装入 CSF、isodata、径向网格及波函数，核对混合记录后调用超精细矩阵例程并关闭报告。 |
| [`rhfs_mpi_main.f90`](rhfs_mpi_main.f90) | rank 0 选择状态与混合来源并广播。集体加载 CSF，准备径向数据和混合向量，运行分布式矩阵元计算并归约各阶段最长耗时。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`hyperfine_one_electron_element.f90`](hyperfine_one_electron_element.f90) | 从轨道 $\kappa,l,j$ 及多极阶生成角选择、旁支 $\kappa$ 和相位，求单电子超精细约化系数。对相关情形同时返回普通和自旋磁矩的角因子。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`hyperfine_radial_integral.f90`](hyperfine_radial_integral.f90) | 在两轨道共同的有效区间内构造 $r^p(P_aQ_b+Q_aP_b)$，乘网格 Jacobian，调用共享径向求积返回标量积分。 |

### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`calculate_hyperfine_matrix_elements.f90`](calculate_hyperfine_matrix_elements.f90) | 按轨道对预制核磁偶极、电四极和 Landé 所需径向、角向因子。枚举合乎 $J$、占据和张量选择的 CSF 对，乘入混合向量形成 ASF 矩阵，再由核矩和角动量因子得到常数及非对角元。 |
| [`calculate_hyperfine_matrix_elements_mpi.f90`](calculate_hyperfine_matrix_elements_mpi.f90) | 以进程数为间隔分配 CSF 列，预筛占据差并按张量阶求矩阵元。将超精细、Landé 及修正数组归约到 rank 0 后进行物理归一化和输出。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`open_summary_file.f90`](open_summary_file.f90) | 根据状态名和混合来源选择 `.h` 或 `.ch`，连接超精细报告单元并检查打开状态。 |
| [`write_hyperfine_energies.f90`](write_hyperfine_energies.f90) | 遍历选定级对，利用能量及宇称标签写允许的非对角偶极／四极结果，并输出供对照的级能量与单位。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 在超精细报告头记录电子、相对论 CSF 和轨道数及相关核数据。 |
