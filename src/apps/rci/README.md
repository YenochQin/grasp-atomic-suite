# rci：相对论组态相互作用

`rci` 在固定的径向轨道上按对称性块组装相对论哈密顿量，求出选定能级的能量和混合系数。
可选修正包括横向 Breit 相互作用、质量位移、真空极化及自能估计。

程序读取 `isodata`、`name_desc.parquet`、`name_header.toml` 和轨道优化得到的 `name.wfn`，
写出 `name.cimix`、`name.csum`、`name.clog` 和续算文件 `rci.res`。
`name.cimix` 可用于生成能级表，也可作为原子性质程序的 CI 输入。

[能级表程序](../rlevels/README.md) · [MPI 版本](../rci_mpi/README.md)

## 计算流程

1. 装入固定径向轨道与描述符基组。根据用户选择建立 Dirac–Coulomb 积分表与可选的 Breit、质量位移、真空极化和自能数据。[程序入口](rci_main.f90) · [计算选项](read_ci_options.f90)
2. 对每个对称性块，逐个 CSF 对把角系数乘到相应的一体、二体径向积分上，再累加启用的修正，按阈值保存下三角条目。已有有效的 `.res` 可以用于恢复已完成的行。[矩阵组装](assemble_hamiltonian_block.f90) · [续算位置恢复](assemble_ci_matrix.f90)
3. 以块平均能量为位移，按矩阵规模和存储选择 `DSPEVX` 或 Davidson 路径。写出选择的根、对称性、块偏移及本征向量。绝对能量由报告中的位移与对应本征值相加。[对角化](diagonalize_hamiltonian_block.f90) · [分块求解与输出](solve_hamiltonian_blocks.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_ci_debug_output.f90`](configure_ci_debug_output.f90) | 为 CI 初始化诊断标志，包括角向输出和调试路径的开关。主流程以这些标志控制后续记录。 |
| [`connect_ci_restart_file.f90`](connect_ci_restart_file.f90) | 为新算例创建 `.res` 并写 `RCI92` 标记，或打开既有续算文件并校验标记。分别调用计算选项或续算读入路径。MPI 按 rank 选择文件。 |
| [`load_csf_basis_data.f90`](load_csf_basis_data.f90) | 调用 `open_csf_blocks` 安装描述符头及块标签，准备按对称性块装入 CSF。矩阵求解循环访问各块。 |
| [`rci_main.f90`](rci_main.f90) | 读取运行模式与状态，打开基组、续算、混合及摘要文件，制备 Slater 与可选 Breit 积分表，再按对称性块求本征对。 |
| [`read_ci_mixing_coefficients.f90`](read_ci_mixing_coefficients.f90) | 按块读取目标根编号，计算每块请求根数、最大根号、根到块的映射和向量总大小，写 `G92MIX` 头供结果输出。 |
| [`read_ci_options.f90`](read_ci_options.f90) | 读取核质量和电荷以及交互的一阶截止、Breit、质量位移、真空极化和自能选项。生成网格与核势并读取径向波函数，再把计算设置写入续算文件头。 |
| [`read_ci_restart.f90`](read_ci_restart.f90) | 恢复核模型、物理修正、径向网格与轨道记录，对照当前基组的电子、轨道和分块计数，并重建径向差分权重。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`accumulate_closed_shell_breit.f90`](accumulate_closed_shell_breit.f90) | 对闭壳层轨道对枚举由角选择允许的张量／径向阶，通过共享 Breit 权重公式形成直接和交换系数，再以带符号的类型码登记闭核贡献。 |
| [`assemble_ci_matrix.f90`](assemble_ci_matrix.f90) | 先检查块大小、变分截止和 rank 布局，扫描已写入的矩阵行并恢复对角和、非零计数与闭核 Breit 位移。从第一条缺行开始继续调用矩阵组装器。 |
| [`assemble_ci_matrix_from_restart.f90`](assemble_ci_matrix_from_restart.f90) | 根据已组装或恢复的统计值求块平均能量和下三角密度。平均值为对角和除以 CSF 数再加闭核 Breit 位移，新算例还输出各积分统计。 |
| [`assemble_hamiltonian_block.f90`](assemble_hamiltonian_block.f90) | 遍历负责的下三角 CSF 行。把 Dirac–Coulomb 角系数与一体、Slater 径向积分相乘，按选项叠加 Breit、质量位移、真空极化和轨道自能，并筛选小非对角元。逐行写稀疏值、列号及闭核能量偏移。 |
| [`collect_breit_angular_coefficient.f90`](collect_breit_angular_coefficient.f90) | 应用 Breit 系数阈值后按需扩容，将左、右轨道四元组、径向阶、积分类型和系数存入矩阵组装缓冲区。 |
| [`diagonalize_hamiltonian_block.f90`](diagonalize_hamiltonian_block.f90) | 读取块的对角和稀疏记录并减去平均能量。根据矩阵大小与存储，选择 LAPACK `DSPEVX` 的打包直接求解，或以主子空间初猜和矩阵乘回调驱动 Davidson 迭代，再保存所选根。 |
| [`initialize_dense_eigenvectors.f90`](initialize_dense_eigenvectors.f90) | 从已打包矩阵提取小主子矩阵，调用 `DSPEVX` 求低端根，按原 CSF 下标回填为全维 Davidson 初猜列。 |
| [`initialize_sparse_eigenvectors.f90`](initialize_sparse_eigenvectors.f90) | 从稀疏行中收集限定子空间内的条目，组装并求解打包主子矩阵，然后把本征向量嵌回大块工作区。 |
| [`multiply_restart_hamiltonian.f90`](multiply_restart_hamiltonian.f90) | 逐批从 `.res` 读取稀疏元素并对输入试向量执行点积与散加，恢复对称矩阵的两侧贡献，按求解模式处理对角项。 |
| [`multiply_sparse_hamiltonian.f90`](multiply_sparse_hamiltonian.f90) | 按列尾偏移遍历内存稀疏结构，在一侧点积的同时将另一侧贡献散加到向量，并加入独立的对角元。 |
| [`satisfies_breit_group1_triangle.f90`](satisfies_breit_group1_triangle.f90) | 依据四个轨道的 $l,j$ 检验第 1 组横向积分允许的三角及宇称组合，作为制表筛选。 |
| [`satisfies_breit_group2_triangle.f90`](satisfies_breit_group2_triangle.f90) | 对第二组横向 Breit 算符使用对应的 $l,j$ 和径向阶三角选择，过滤禁止预制的轨道组合。 |
| [`satisfies_coulomb_triangle.f90`](satisfies_coulomb_triangle.f90) | 对构成 Slater 项的两个轨道对执行角动量三角和宇称条件检查，决定是否枚举该多极阶。 |
| [`solve_hamiltonian_blocks.f90`](solve_hamiltonian_blocks.f90) | 遍历各块，加载 CSF、设置目标根并分配工作区。根据可用续算记录恢复或组装下三角矩阵，求所选本征对，将块平均能量、根和向量写入混合文件，并输出可选自能结果。 |

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
| [`electron_screening_factor.f90`](electron_screening_factor.f90) | 按轨道能量、$n$ 和 $\kappa$ 反求等效 Coulomb 电荷，返回 $Z-Z_{\mathrm{eff}}$。遇到无实解的试探值时返回越界标记。 |
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
| [`tabulate_breit_integrals_2.f90`](tabulate_breit_integrals_2.f90) | 针对第二类 Breit 角选择生成各阶的轨道四元组，调用径向内核并保存到第二组积分值和索引表。 |
| [`tabulate_self_energy_integrals.f90`](tabulate_self_energy_integrals.f90) | 对截止 $n$ 内的轨道，用近核密度比乘氢样自能因子及 $Z^4/(\pi n^3c^3)$ 制表。超过主量子数上限的轨道自能项置零。 |
| [`tabulate_slater_integrals.f90`](tabulate_slater_integrals.f90) | 按多极阶枚举轨道四元组，用三角与对称条件筛选。先确定每阶的表范围，再以规范打包键保存 Coulomb 径向积分。 |
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
| [`open_mixing_output.f90`](open_mixing_output.f90) | 打开该状态的 CI 混合文件，进入目标根选择以及 `G92MIX` 头和维数记录的写出。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 打开状态的 `.csum` 文件并设置共享摘要单元，供选项、矩阵规模和能级结果写出。 |
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
