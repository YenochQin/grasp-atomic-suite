# rci_mpi：并行相对论组态相互作用

`rci_mpi` 在固定的径向轨道上求解相对论组态相互作用，将哈密顿量组装和稀疏矩阵与向量的乘积
分配给多个 MPI 进程。

新算例读取 `isodata`、`name_desc.parquet`、`name_header.toml` 和 `name.wfn`。
启用零阶优先分区时从头文件读取 `[block_info].zero_order_counts`。
程序在永久 I/O 目录写出混合系数 `name.cimix` 和摘要 `name.cisum`。rank 0 在自己的工作目录写出日志
`name.cilog`，各 rank 在各自工作目录写出 `name.rci.rankNNN.restart` 续算文件。

[串行版本与可选修正](../rci/README.md) · [原生 CSF 输入及零阶优先分区](../../../docs/implemented/design-0001-native-csf-descriptors.md)

## 计算流程

1. MPI 初始化后，集体打开描述符块和各 rank 的本地 `.rci.rankNNN.restart`。rank 0 读取目标根和物理选项。零阶优先选项从描述符头取各块截止位置。[入口](rci_mpi_main.f90) · [分区选项](read_ci_options.f90)
2. 各 rank 分担预制的 Slater 与 Breit 积分，再合并径向表。分块加载 CSF 后，各 rank 组装负责的稀疏行，并归约对角和以确定共同的能量位移。[Slater 表归约](tabulate_slater_integrals_mpi.f90) · [块组装](assemble_hamiltonian_block.f90) · [能量位移](assemble_ci_matrix_from_restart.f90)
3. 在迭代求解路径中，每个 rank 计算本地的矩阵乘积，再归约得到 Davidson 所用的完整向量。所选能级及本征向量最终由 rank 0 写入共享混合文件。[矩阵乘子程序](multiply_sparse_hamiltonian.f90) · [求根及写出](solve_hamiltonian_blocks.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_ci_debug_output.f90`](configure_ci_debug_output.f90) | 为 CI 初始化诊断标志，包括角向输出和调试路径的开关。主流程以这些标志控制后续记录。 |
| [`connect_ci_restart_file.f90`](connect_ci_restart_file.f90) | 按 rank 打开 `.rci.rankNNN.restart`。每个 rank 在续算时校验自己的 `R92RS2` 分片及其 V3 描述符 SHA-256。旧 `R92RES` 分片须通过冷启动重新生成。 |
| [`rci_mpi_main.f90`](rci_mpi_main.f90) | 初始化各 rank 的永久输入和临时目录，广播设置、打开分块描述符及本地续算流。并行生成径向积分表后进入分块求解，由 rank 0 写共享报告。 |
| [`read_ci_mixing_coefficients_mpi.f90`](read_ci_mixing_coefficients_mpi.f90) | rank 0 交互取得目标根编号，广播根数、最大根号及选择列表。各 rank 由相同块维数计算本征向量布局，rank 0 写混合头。 |
| [`read_ci_options.f90`](read_ci_options.f90) | rank 0 读取非默认速度以及各物理修正选项并广播。若启用零阶优先分区，从描述符头取得每块截止。读取核和径向数据后写本地续算设置。 |
| [`read_ci_restart.f90`](read_ci_restart.f90) | 恢复核模型、物理修正、径向网格与轨道记录，对照当前基组的电子、轨道和分块计数，并重建径向差分权重。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`accumulate_closed_shell_breit.f90`](accumulate_closed_shell_breit.f90) | 对闭壳层轨道对枚举由角选择允许的张量／径向阶，通过共享 Breit 权重公式形成直接和交换系数，再以带符号的类型码登记闭核贡献。 |
| [`assemble_ci_matrix.f90`](assemble_ci_matrix.f90) | 先检查块大小、变分截止和 rank 布局，扫描已写入的矩阵行并恢复对角和、非零计数与闭核 Breit 位移。从第一条缺行开始继续调用矩阵组装器。 |
| [`assemble_ci_matrix_from_restart.f90`](assemble_ci_matrix_from_restart.f90) | 在新矩阵路径把积分计数归约给 rank 0，对每 rank 存储项数取最大值。对角能量求和后除以全块 CSF 数，加上广播的闭核 Breit 位移作为块能量偏移。 |
| [`assemble_hamiltonian_block.f90`](assemble_hamiltonian_block.f90) | 遍历负责的下三角 CSF 行。把 Dirac–Coulomb 角系数与一体、Slater 径向积分相乘，按选项叠加 Breit、质量位移、真空极化和轨道自能，并筛选小非对角元。逐行写稀疏值、列号及闭核能量偏移。 |
| [`collect_breit_angular_coefficient.f90`](collect_breit_angular_coefficient.f90) | 应用 Breit 系数阈值后按需扩容，将左、右轨道四元组、径向阶、积分类型和系数存入矩阵组装缓冲区。 |
| [`diagonalize_hamiltonian_block.f90`](diagonalize_hamiltonian_block.f90) | 读取块的对角和稀疏记录并减去平均能量。根据矩阵大小与存储，选择 LAPACK `DSPEVX` 的打包直接求解，或以主子空间初猜和矩阵乘回调驱动 Davidson 迭代，再保存所选根。 |
| [`initialize_dense_eigenvectors.f90`](initialize_dense_eigenvectors.f90) | 从已打包矩阵提取小主子矩阵，调用 `DSPEVX` 求低端根，按原 CSF 下标回填为全维 Davidson 初猜列。 |
| [`initialize_sparse_eigenvectors.f90`](initialize_sparse_eigenvectors.f90) | 从稀疏行中收集限定子空间内的条目，组装并求解打包主子矩阵，然后把本征向量嵌回大块工作区。 |
| [`multiply_restart_hamiltonian.f90`](multiply_restart_hamiltonian.f90) | 逐批从 `.rci.rankNNN.restart` 读取稀疏元素并对输入试向量执行点积与散加，恢复对称矩阵的两侧贡献，按求解模式处理对角项。 |
| [`multiply_sparse_hamiltonian.f90`](multiply_sparse_hamiltonian.f90) | 在本 rank 的稀疏列上执行对称散加与点积并处理本地对角元，随后对全部试向量的乘积作 MPI 求和。 |
| [`satisfies_breit_group1_triangle.f90`](satisfies_breit_group1_triangle.f90) | 依据四个轨道的 $l,j$ 检验第 1 组横向积分允许的三角及宇称组合，作为制表筛选。 |
| [`satisfies_breit_group2_triangle.f90`](satisfies_breit_group2_triangle.f90) | 对第二组横向 Breit 算符使用对应的 $l,j$ 和径向阶三角选择，过滤禁止预制的轨道组合。 |
| [`satisfies_coulomb_triangle.f90`](satisfies_coulomb_triangle.f90) | 对构成 Slater 项的两个轨道对执行角动量三角和宇称条件检查，决定是否枚举该多极阶。 |
| [`solve_hamiltonian_blocks.f90`](solve_hamiltonian_blocks.f90) | 各 rank 集体取当前 CSF 块，生成各自的矩阵行并汇总块能量，进入分布式矩阵乘或打包直接对角化路径。rank 0 写混合、能级及选配自能结果。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`breit_radial_integral.f90`](breit_radial_integral.f90) | 按 Breit 类型确定两侧选用的大、小分量、交换轨道顺序和阶次，调用对应的横向径向求积构成四轨道积分。 |
| [`breit_radial_integral_1.f90`](breit_radial_integral_1.f90) | 按轨道四元组打包键，在 `first_breit_integral_group_offsets` 给出的阶次区间做二分查找，返回预制的第 1 类积分。 |
| [`breit_radial_integral_2.f90`](breit_radial_integral_2.f90) | 在第 2 类 Breit 表的指定多极阶区间二分查找轨道四元组，检查索引匹配后返回已制表积分。 |
| [`breit_radial_integral_3.f90`](breit_radial_integral_3.f90) | 以径向阶和四个轨道打包第 3 类键。命中读取第 3 类表，缺项调用对应 Breit 内核并按序插入。 |
| [`breit_radial_integral_4.f90`](breit_radial_integral_4.f90) | 查找包含径向阶的第 4 类轨道键，未命中时计算第 4 类径向积分并保存到有序缓存。 |
| [`breit_radial_integral_5.f90`](breit_radial_integral_5.f90) | 对第 5 类实现独立的有序查找与按需计算。插入记录时同步移动索引和值。 |
| [`breit_radial_integral_6.f90`](breit_radial_integral_6.f90) | 打包第 6 类径向阶和轨道键，二分查表，缺项计算后在第 6 类的索引与数值数组中有序保存。 |
| [`build_breit_bessel_kernels.f90`](build_breit_bessel_kernels.f90) | 对频率缩放后的径向参数，按量级用级数或递推生成球 Bessel／Neumann 核，并保存供横向积分复用的缩放值。 |
| [`build_fourth_order_vacuum_polarization.f90`](build_fourth_order_vacuum_polarization.f90) | 计算点核四阶势，对需处理的有限核部分用 Fullerton–Rinker $L_n$ 核和核电荷密度做径向积分，保存四阶势。 |
| [`build_nuclear_charge_distribution.f90`](build_nuclear_charge_distribution.f90) | 由半密度半径和扩散参数形成 Fermi 核密度，用交替指数级数求归一化并在径向网格上制表。 |
| [`build_second_order_vacuum_polarization.f90`](build_second_order_vacuum_polarization.f90) | 用 Fullerton–Rinker $K_n$ 核构造点核表达式。对有限核用径向求积卷积核密度，返回二阶势表。 |
| [`build_vacuum_polarization.f90`](build_vacuum_polarization.f90) | 设置当前核分布数据，依次构造二阶和四阶真空极化势，供一体真空极化矩阵元合用。 |
| [`calculate_vacuum_polarization_integral.f90`](calculate_vacuum_polarization_integral.f90) | 对规范化的一体轨道对查真空极化积分缓存。缺项时调用径向势积分，再按序保存键与值。 |
| [`calculate_vinti_integral.f90`](calculate_vinti_integral.f90) | 将 Vinti 轨道对映射到保留方向的键，二分查有序表，首次访问时计算并缓存对应的 SMS 径向积分。 |
| [`coulomb_radial_integral.f90`](coulomb_radial_integral.f90) | 构造一对轨道的径向乘积，用 Hartree 内积分生成辅助函数，再与第二对轨道径向乘积及 Jacobian 相乘，求积形成四轨道 Coulomb 积分。 |
| [`evaluate_breit_hartree_z.f90`](evaluate_breit_hartree_z.f90) | 逐区间积分带 $(s/r)^k$ 权重的径向源，累积从原点至当前点的 Hartree $Z$ 函数，写回供外积分使用的径向工作数组。 |
| [`fullerton_rinker_k_function.f90`](fullerton_rinker_k_function.f90) | 按输入 $x$ 的区间选择小参数解析表达或较大参数的逼近，返回二阶真空极化求积需要的 $K_n(x)$。 |
| [`fullerton_rinker_l_function.f90`](fullerton_rinker_l_function.f90) | 在不同 $x$ 区间使用 Fullerton–Rinker 的解析级数或衰减逼近，返回四阶势需要的 $L_n(x)$。 |
| [`hydrogenic_radial_overlap.f90`](hydrogenic_radial_overlap.f90) | 用给定有效电荷构造同 $n,\kappa$ 的解析氢样轨道，在有效网格上积分 $P P_H+Q Q_H$。 |
| [`hydrogenic_self_energy_factor.f90`](hydrogenic_self_energy_factor.f90) | 低主量子数用 Mohr 表。较高 $n$ 的指定 $\kappa$ 借用 $n=2$ 的 Mohr 项，其他轨道进入 Klarsfeld–Maquet 近似。 |
| [`kinetic_energy_integral.f90`](kinetic_energy_integral.f90) | 查找当前轨道对的质量位移动能积分。未命中时求当前径向轨道的动能矩阵元并插入有序表。 |
| [`klarsfeld_maquet_self_energy_factor.f90`](klarsfeld_maquet_self_energy_factor.f90) | 由 $n,\kappa$ 查 Bethe 对数并求 $Z\alpha$ 下的近似展开项，组合成氢样自能因子。 |
| [`lookup_breit_radial_integral.f90`](lookup_breit_radial_integral.f90) | 按积分类型与多极阶准备对应的 Bessel／Neumann 核，以轨道对和相邻阶次生成核缓存键，再将径向分量组合传给横向积分内核。 |
| [`lookup_slater_integral.f90`](lookup_slater_integral.f90) | 将两轨道对归一到规范次序，用多极阶偏移限定二分区间，从预制 Slater 表中读取对应径向积分。 |
| [`mohr_self_energy_factor.f90`](mohr_self_energy_factor.f90) | 按 $n$ 和 $\kappa$ 选取 Mohr 表的序列，以核电荷为自变量作 Aitken 插值得到无量纲因子。 |
| [`one_electron_integral.f90`](one_electron_integral.f90) | 把一体轨道对归一化并打包为有序键，命中返回缓存。缺项时计算 Dirac–Coulomb 一体积分并插入，必要时扩充表容量。 |
| [`prepare_correction_integral_tables.f90`](prepare_correction_integral_tables.f90) | 初始化各类积分的计数与首次查找标志，由核质量得到反质量因子。启用真空极化时先制备核密度和相应势数组。 |
| [`radial_density_ratio.f90`](radial_density_ratio.f90) | 为同一 $n,\kappa$ 计算氢样参考轨道，在选定近核积分范围内分别累计数值及氢样轨道的 $P^2+Q^2$，返回两积分之比供自能缩放。 |
| [`tabulate_breit_integrals_1.f90`](tabulate_breit_integrals_1.f90) | 在第一类 Breit 角选择下按多极阶枚举轨道组合，建立各阶偏移并填入打包键和第 1 类径向积分值。 |
| [`tabulate_breit_integrals_1_mpi.f90`](tabulate_breit_integrals_1_mpi.f90) | 分担第 1 类 Breit 积分表计算，MPI 全归约合并分派的打包键与实数结果。 |
| [`tabulate_breit_integrals_2.f90`](tabulate_breit_integrals_2.f90) | 针对第二类 Breit 角选择生成各阶的轨道四元组，调用径向内核并保存到第二组积分值和索引表。 |
| [`tabulate_breit_integrals_2_mpi.f90`](tabulate_breit_integrals_2_mpi.f90) | 分担第 2 类 Breit 积分表计算，MPI 全归约合并分派的打包键与实数结果。 |
| [`tabulate_self_energy_integrals.f90`](tabulate_self_energy_integrals.f90) | 对截止 $n$ 内的轨道，用近核密度比乘氢样自能因子及 $Z^4/(\pi n^3c^3)$ 制表。超过主量子数上限的轨道自能项置零。 |
| [`tabulate_slater_integrals.f90`](tabulate_slater_integrals.f90) | 按多极阶枚举轨道四元组，用三角与对称条件筛选。先确定每阶的表范围，再以规范打包键保存 Coulomb 径向积分。 |
| [`tabulate_slater_integrals_mpi.f90`](tabulate_slater_integrals_mpi.f90) | 调用 Slater 制表例程计算分派的积分值，用 MPI 整数和实数求和合并键及数值数组。 |
| [`transverse_radial_integral.f90`](transverse_radial_integral.f90) | 对指定的轨道乘积应用频率相关的径向核与 Hartree 内积分，再在有限的重叠网格上完成外积分。内部按参数选择允许的近似路径。 |
| [`vacuum_polarization_integral.f90`](vacuum_polarization_integral.f90) | 将两轨道的大、小分量乘积与二阶、四阶势相乘并在共同有效径向区间求积，得到一体真空极化矩阵元。 |
| [`vinti_radial_integral.f90`](vinti_radial_integral.f90) | 差分求径向导数，按两轨道 $\kappa$ 对应的轨道角动量构造 $1/r$ 项，将大、小分量合成 Vinti 被积函数并求积。 |

### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`build_qed_self_energy.f90`](build_qed_self_energy.f90) | 对某个本征态按 CSF 占据乘 $c_i^2$ 得到各轨道的广义占据。同时用实际核电荷及近核密度比构造供态自能估计的轨道修正。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`open_mixing_output_mpi.f90`](open_mixing_output_mpi.f90) | 由 rank 0 打开 CI 混合文件。各 rank 共同进入能级选择广播，文件头仅由输出 rank 写入。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 打开状态的 `.cisum` 文件并设置共享摘要单元，供选项、矩阵规模和能级结果写出。 |
| [`read_ci_level_selection.f90`](read_ci_level_selection.f90) | 扫描已求出的每个态的 CSF 系数，用有序插入保留绝对值最大的至多五项，然后将带符号系数与块内 CSF 索引写入总结。 |
| [`write_ci_eigenvectors.f90`](write_ci_eigenvectors.f90) | 把所选根和对应本征向量按每页八列展开输出，每行对应一个 CSF，供诊断完整混合系数。 |
| [`write_hamiltonian_records.f90`](write_hamiltonian_records.f90) | 逐条读续算流里的矩阵行，解析非零值和行列索引，将 $H_{ij}$ 写到调试单元供核对矩阵。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 写入 CSF 分块、零阶／一阶截断、核模型、径向网格及启用的哈密顿量修正，记录各块的目标根与轨道参数。 |
| [`write_level_energies.f90`](write_level_energies.f90) | 从块偏移与本征值计算各级总能量，按输出模式追加相邻能级间隔和相对块最低根的能量，并换算单位。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`ci_restart_units.f90`](ci_restart_units.f90) | 保存当前打开的矩阵续算逻辑单元，供组装、稀疏乘法和对角化按同一文件定位。 |
| [`find_breit_integral_index.f90`](find_breit_integral_index.f90) | 根据积分类型选择六个 Breit 索引数组之一，按给定表位置返回该条目的打包标签。 |
| [`reduce_sparse_matrix_product.f90`](reduce_sparse_matrix_product.f90) | 计算一段非对角项的点积 `dot_product(da,db)`，并原地更新另一侧向量 `dc += dconst*da`，为对称矩阵乘提供局部操作。 |
