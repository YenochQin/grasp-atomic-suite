# rmcdhf_mpi：并行 MCDHF

`rmcdhf_mpi` 将 MCDHF 哈密顿量和轨道方程的计算分配给多个 MPI 进程。
每个 rank 在迭代中从磁盘读取自己的 MCP 通道。

使用与 `rangular_mpi` 相同的基组、进程数和临时目录布局，并在永久 I/O 目录提供
`isodata`、状态的描述符文件对及 `name.initial.wfn`。程序将 `name.wfn`、SCF 摘要和
OL/EOL 计算的 `name.mix` 写入永久 I/O 目录，SCF 日志写入启动目录。

需要通过 `npol.inp` 启用核极化时，可选用 MPI 内存版。
[角系数程序](../rangular_mpi/README.md) · [MPI 内存版](../rmcdhf_mem_mpi/README.md) · [MPI 诊断](../../../docs/implemented/design-0021-rmcdhf-mpi-diagnostics.md)

势系数复用当前密度权重下的收缩结果。[格式与诊断选项](../../../docs/implemented/design-0020-rmcdhf-mcp-io.md)

## 计算流程

1. 同步状态与网格设置，各 rank 读取并核对自己的 MCP 文件。集体加载基组与初猜，建立相同的核势和 SCF 权重数据。[MPI 输入](rmcdhf_mpi_main.f90) · [径向与 SCF 设置](read_scf_options_mpi.f90)
2. 各 rank 对本地稀疏记录计算径向积分和哈密顿量贡献，归约块对角和以取得共同能量偏移。Davidson 矩阵乘根据本地稀疏矩阵计算并合并各 rank 的乘积。[分块哈密顿量](solve_hamiltonian_blocks_mpi.f90) · [分布式求根](diagonalize_hamiltonian_block_mpi.f90)
3. 归约广义占据和径向势，计算同 $\kappa$ 轨道的拉格朗日乘子，执行一致的轨道改善、阻尼和正交化。两个轨道均可变的 MPI 乘子分支使用势和公式。[乘子](calculate_lagrange_multipliers_mpi.f90) · [轨道改善](improve_orbital_mpi.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 初始化 SCF 径向、角向、矩阵与正交化诊断位，按需打开指定 `.dbg` 路径。 |
| [`configure_debug_output_mpi.f90`](configure_debug_output_mpi.f90) | 在 rank 0 处理调试选项并广播诊断位，使各 rank 采用相同的临时诊断设置。 |
| [`connect_mcp_files.f90`](connect_mcp_files.f90) | 打开本 rank 的 `mcp.mcpidx`，校验分片身份及布局并安装基组和相互作用元数据。 |
| [`load_csf_basis_data_mpi.f90`](load_csf_basis_data_mpi.f90) | 调用 MPI 适配器，由 rank 0 验证并读取描述符，各 rank 对照自己的 MCP 元数据并取得广播的占据和块信息。 |
| [`read_al_average_options_mpi.f90`](read_al_average_options_mpi.f90) | 由 rank 0 选择 AL/EAL 权重和迭代控制选项并广播，使各 rank 构造一致的平均占据和正交化设置。 |
| [`read_al_level_weights.f90`](read_al_level_weights.f90) | 按等权、$2J+1$ 或用户给出的正权重填入向量，对总和归一化供平均能量与占据计算。 |
| [`read_eol_average_options_mpi.f90`](read_eol_average_options_mpi.f90) | rank 0 读取目标根和轨道优化控制，广播各块最高根、请求根列表、权重与更新顺序，为共享波函数更新分配状态。 |
| [`read_eol_level_weights.f90`](read_eol_level_weights.f90) | 等权与统计权用标记保留到根求出后确定。自定义输入要求正权重，立即除以权重总和。 |
| [`read_scf_options_mpi.f90`](read_scf_options_mpi.f90) | 由 rank 0 读核、初猜与交互选项，广播网格、物理和 SCF 控制数据。各 rank 准备核势与等价的初始径向分量，然后选择 AL/EAL 或 OL/EOL 的目标及权重。 |
| [`rmcdhf_mpi_main.f90`](rmcdhf_mpi_main.f90) | 初始化 rank 工作目录，由 rank 0 读取状态和永久输入路径。各 rank 验证自己的 MCP 记录，集体载入描述符并初始化 SCF，rank 0 保存轨道和相应级的结果。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`assemble_hamiltonian_block.f90`](assemble_hamiltonian_block.f90) | 先按各 CSF 的子壳层占据累加对角直接项。循 MCP 的 T 和各阶 V 标签计算一体／二体径向积分，以 `角系数 × 积分` 散加到记录指定的稀疏矩阵位置。 |
| [`assemble_orbital_potential_coefficients.f90`](assemble_orbital_potential_coefficients.f90) | 读一体和各多极二体 MCP 条目，利用当前 CSF 密度矩阵将贡献收缩。按被优化轨道及积分标签分类，形成直接势、交换源和非对角项的系数列表。 |
| [`csf_density_matrix_element.f90`](csf_density_matrix_element.f90) | 用块偏移把局部 CSF 序号映射到混合数组。EOL 返回 $D_{ij}=\sum_s w_sc_{is}c_{js}$。平均层路径仅返回加权对角元。 |
| [`diagonalize_hamiltonian_block_mpi.f90`](diagonalize_hamiltonian_block_mpi.f90) | 从本地矩阵各取对角并归约，构造分布式 Davidson 初猜。回调合并各 rank 的稀疏乘积，返回同一组目标根、位移能量和向量。 |
| [`direct_angular_coefficient.f90`](direct_angular_coefficient.f90) | 零阶直接项按同壳电子对数 $q(q-1)/2$ 或两壳占据乘积计算。同壳的非零阶闭壳项按相对论 3j 系数计算。 |
| [`exchange_angular_coefficient.f90`](exchange_angular_coefficient.f90) | 当至少一个子壳层闭合时，以占据数乘积乘相对论 3j 系数平方并取负，形成对应多极阶的交换系数。 |
| [`run_scf_iterations_mpi.f90`](run_scf_iterations_mpi.f90) | 同步迭代和轨道选择设置，循环执行分布式矩阵组装／求根、权重与乘子更新和轨道改善。用一致的收敛判据结束循环，由输出 rank 写轨道与最终摘要。 |
| [`solve_hamiltonian_blocks_mpi.f90`](solve_hamiltonian_blocks_mpi.f90) | 各 rank 用本地 MCP 组装块矩阵，对对角和做全归约并移去统一的平均值，然后分布式求根。更新与上一轮对应的混合向量，写所选块的能量和系数。 |
| [`update_scf_mixing_mpi.f90`](update_scf_mixing_mpi.f90) | 依据目标根生成和归一化权重，以本 rank 的 CSF 计算轨道占据的局部加和，用 MPI 归约获得完整占据后计算加权态能量。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`apply_lagrange_constraints.f90`](apply_lagrange_constraints.f90) | 扫描与当前轨道有关的正交约束，按占据或代码给定的归一化权重选取乘子，将约束轨道的径向分量加入对应源项。 |
| [`apply_off_diagonal_integrals.f90`](apply_off_diagonal_integrals.f90) | 对非对角一电子积分的来源轨道计算径向导数，以能量、核势、$\kappa$ 和系数构造耦合项，更新大小分量交换源。 |
| [`build_direct_radial_potential.f90`](build_direct_radial_potential.f90) | 以核势为起点，根据直接势系数记录提取每个轨道对，逐阶计算 Hartree $Y^k$ 径向函数并按系数累加到当前轨道的 $Y(r)$ 势。 |
| [`build_exchange_radial_potential.f90`](build_exchange_radial_potential.f90) | 解包每项的径向积分轨道与源轨道，计算 Hartree 相互作用，将它与源轨道的 $P,Q$ 乘积分别累加进大小分量交换源。 |
| [`build_orbital_potential_mpi.f90`](build_orbital_potential_mpi.f90) | 各 rank 按本地角系数生成直接、交换及非对角径向贡献，再将分布式势与源项归约为轨道求解使用的完整数组。 |
| [`build_radial_inhomogeneous_terms.f90`](build_radial_inhomogeneous_terms.f90) | 以交换源、网格因子与步长构造 $P,Q$ 的非齐次右端并形成相邻步进需要的源数组。 |
| [`build_variation_source_term.f90`](build_variation_source_term.f90) | 以 $-Q$ 与 $P$ 分别乘半步长、Jacobian 和 $1/c$，累加相邻网格的右端项，形成随轨道能量变化的线性响应源。 |
| [`calculate_deferred_corrections.f90`](calculate_deferred_corrections.f90) | 从第三个活动网格点开始，对上一轮 $P,Q$ 施加权重为 `1,-15,40` 除以 `120` 的差分组合，得到径向步进的延迟修正，边界处置零。 |
| [`calculate_lagrange_multipliers_mpi.f90`](calculate_lagrange_multipliers_mpi.f90) | 分布式生成势和一电子修正，归约约束所需的积分。两个轨道均可变的分支使用势和与占据倒数和的公式，再记录乘子。 |
| [`check_orbital_damping.f90`](check_orbital_damping.f90) | 比较连续两次相对能量变化。启用自适应阻尼时，若乘积小于 `-1e-4` 则加强阻尼，否则将阻尼减半，并保存本次变化。 |
| [`clear_radial_inhomogeneous_terms.f90`](clear_radial_inhomogeneous_terms.f90) | 把当前大小分量的两个径向源数组置零，使下一次积分先按齐次方程推进。 |
| [`damp_orbital_update.f90`](damp_orbital_update.f90) | 按阻尼因子混合新旧大、小分量与原点参数，处理有效网格长度和符号，重新归一化并检查径向节点。 |
| [`energy_outside_bounds.f90`](energy_outside_bounds.f90) | 将候选能量与当前上下界做严格开区间比较，落在区间外或等于端点时通知更新例程采取修正。 |
| [`estimate_orbital_energy.f90`](estimate_orbital_energy.f90) | 由 $n,\kappa$ 与选取的核电荷构造氢样 Dirac 能量尺度，建立当前轨道允许的试探能量上下界。 |
| [`evaluate_radial_self_consistency.f90`](evaluate_radial_self_consistency.f90) | 在新旧轨道共有的有效网格上取 $\max_r(\lvert\Delta P\rvert+\lvert\Delta Q\rvert)$，再乘广义占据的平方根，写入当前轨道的自洽误差。 |
| [`find_least_consistent_orbital.f90`](find_least_consistent_orbital.f90) | 在可变轨道列表中查找最大的自洽误差，返回轨道索引及误差供自适应更新顺序和停止条件使用。 |
| [`improve_orbital_mpi.f90`](improve_orbital_mpi.f90) | 用归约后的径向势执行轨道求解，协调各 rank 的系数缓冲和更新控制。对试探解做归一化、阻尼及正交化，并记录共同误差。 |
| [`integrate_radial_inward.f90`](integrate_radial_inward.f90) | 以前向消元建立尾部五带系统的三角因子，按尾振幅与网格边界选取外端，再逆向回代至接点。 |
| [`integrate_radial_outward.f90`](integrate_radial_outward.f90) | 从已给出的前三个原点样本开始，用变换网格的步长、$\kappa r'/r$ 及邻点源项逐步推进耦合 $P,Q$，一直到内外解接点。 |
| [`orthogonalize_same_kappa_orbitals.f90`](orthogonalize_same_kappa_orbitals.f90) | 按固定／可变轨道策略从同 $\kappa$ 轨道减去相互投影，随径向分量一起变换原点系数，再做归一化及有效网格处理。 |
| [`orthogonalize_to_more_consistent_orbitals.f90`](orthogonalize_to_more_consistent_orbitals.f90) | 选择同 $\kappa$ 且优化误差更小的参考轨道，计算重叠并减去投影。检查归一化及节点计数，才使用经修正的试探解。 |
| [`solve_radial_eigenvalue.f90`](solve_radial_eigenvalue.f90) | 用当前数值轨道分别积分交换、直接势、小分量静质量、$\kappa/r$ 和导数项，将五项相加得到轨道能量参数。 |
| [`solve_radial_orbital.f90`](solve_radial_orbital.f90) | 为试探能量建立径向系数和原点条件，先求齐次解再按需加入交换／约束源。分别由内外积分到接点，检查节点数、能量括区和接点失配，用变分源修正能量并归一化试探轨道。 |
| [`update_trial_orbital_energy.f90`](update_trial_orbital_energy.f90) | 根据当前节点数、目标节点数和能量修正移动括区。选择限制步长、区间内插或重新估计能量，并传播找不到允许能量时的状态。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`finish_summary.f90`](finish_summary.f90) | 输出最后一轮的能量、能级权重和组态成分，按选取的平均模式完成报告并关闭相应输出单元。 |
| [`open_mixing_output.f90`](open_mixing_output.f90) | 为 OL/EOL 创建混合文件并写 `G92MIX` 标记和预备维数，供每轮矩阵求解写各块系数。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 将状态摘要路径连接到 SCF 文本单元，检查是否成功打开，供初始设置和结束时的能量表写入。 |
| [`print_radial_wavefunction.f90`](print_radial_wavefunction.f90) | 按列打印网格上的大、小分量及相关势的诊断样本，可调用 ASCII 绘图例程查看径向形状。 |
| [`print_scf_results.f90`](print_scf_results.f90) | 把轨道编号转换为谱学标签，列出将优化的轨道及其控制选项，供交互核对轨道选择。 |
| [`write_csf_weights.f90`](write_csf_weights.f90) | 据每块 CSF、根和混合系数偏移扫过本征向量，选出幅值最大的组态分量，以块内 CSF 序号、符号和对称性输出。 |
| [`write_hamiltonian_records_mpi.f90`](write_hamiltonian_records_mpi.f90) | 遍历当前 rank 负责的稀疏列和对角索引，将所选的矩阵元输出到调试单元。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 记录核模型、离化电荷、基组规模、径向网格、优化轨道和目标能级设置，建立计算报告的输入部分。 |
| [`write_scf_energies.f90`](write_scf_energies.f90) | 对请求的根恢复块平均能量偏移，按 `2J+1` 和宇称标注级能量，打印 Hartree、cm⁻¹ 和 eV 及当前模式要求的能差。 |
| [`write_scf_radial_wavefunctions.f90`](write_scf_radial_wavefunctions.f90) | 按 `GRWFRAD1` 写优化轨道和共享径向网格。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`basis_csf_parity.f90`](basis_csf_parity.f90) | 把 CSF 序号交给 RMCDHF 基组查询接口，由占据／块标签状态读取该 CSF 宇称。 |
| [`basis_csf_two_j_plus_one.f90`](basis_csf_two_j_plus_one.f90) | 从 CSF 对应的块状态取得角动量编码 $2J+1$，供统计权重和能级输出使用。 |
