# rmcdhf_mem：内存版 MCDHF

`rmcdhf_mem` 串行求解 MCDHF，并在内存中缓存 MCP 角系数及其索引。
自洽场迭代组装哈密顿量和轨道方程时复用这些数据。
缓存装载时，若主机提供 Linux 的 `/proc/meminfo`，日志还会打印系统总内存与空闲内存。缓存分配量由 MCP 记录量决定。

程序读取 `isodata`、状态的描述符文件对、`name.initial.wfn` 和 `rangular` 生成的
`mcp.NN` 通道，支持通过 `npol.inp` 启用核极化。输出包括 `name.wfn`、SCF 摘要和日志，
以及 OL/EOL 计算的 `name.mix`。

[串行版本与状态文件](../rmcdhf/README.md) · [MCP 存储及实现差异](../../../docs/reviews/evidence-0007-rmcdhf-memory-variant-comparison.md)

## 计算流程

1. 连接 MCP 文件时先统计每个通道的记录量与分块边界，按实际需求分配空间并装入角系数、打包积分标签及稀疏结构。[MCP 缓存](cache_mcp_files.f90)
2. 每轮从缓存取块结构和角系数，按当前径向轨道代入积分，组装并对角化块哈密顿量。将目标态混合向量按权重收缩成轨道占据和势的系数。[块求解](solve_hamiltonian_blocks.f90) · [势系数收缩](assemble_orbital_potential_coefficients.f90)
3. 构造直接和交换项，求解带约束的径向 Dirac 方程并更新轨道。尾部向内积分以振幅衰减判据选外端点，随后回代到接点。SCF 循环根据轨道误差和加权能量变化停止并写出结果。[径向尾部](integrate_radial_inward.f90) · [SCF 迭代](run_scf_iterations.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 初始化 SCF 径向、角向、矩阵与正交化诊断位，按需打开指定 `.dbg` 路径。 |
| [`connect_mcp_files.f90`](connect_mcp_files.f90) | 读 `mcp.30` 的闭核数、块长度、标签和最大多极阶，再核对各通道的电子、CSF、轨道数及 rank／进程数。保留各块位置与相互作用截止。内存版在这里预载已检查的通道。 |
| [`load_csf_basis_data.f90`](load_csf_basis_data.f90) | 以已经读出的 MCP 元数据调用 RMCDHF 描述符适配器，核查轨道、CSF 和对称性块，随后装入占据及块量子数，供 SCF 用统一的全基组下标访问。 |
| [`read_al_average_options.f90`](read_al_average_options.f90) | 设置 AL/EAL 的权重、迭代上限及正交化策略，将 CSF 占据按选定权重收缩成各轨道的平均占据。 |
| [`read_al_level_weights.f90`](read_al_level_weights.f90) | 按等权、$2J+1$ 或用户给出的正权重填入向量，对总和归一化供平均能量与占据计算。 |
| [`read_eol_average_options.f90`](read_eol_average_options.f90) | 逐块读取需优化的能级并登记块偏移。设置固定和关联轨道、更新次序、能级权重及混合／轨道阻尼。 |
| [`read_eol_level_weights.f90`](read_eol_level_weights.f90) | 等权与统计权用标记保留到根求出后确定。自定义输入要求正权重，立即除以权重总和。 |
| [`read_scf_options.f90`](read_scf_options.f90) | 读取 `isodata`、网格及收敛设置，初始化核 Coulomb 势并应用可选核极化。载入初猜径向轨道，随后设置变分轨道、能级、权重、阻尼和迭代上限。 |
| [`rmcdhf_mem_main.f90`](rmcdhf_mem_main.f90) | 验证状态和 MCP 数据，并在连接各通道时将记录读入缓存。核对描述符、核数据和初猜，准备摘要及混合输出后开始串行轨道优化。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`assemble_hamiltonian_block.f90`](assemble_hamiltonian_block.f90) | 先按各 CSF 的子壳层占据累加对角直接项。循 MCP 的 T 和各阶 V 标签计算一体／二体径向积分，以 `角系数 × 积分` 散加到记录指定的稀疏矩阵位置。 |
| [`assemble_orbital_potential_coefficients.f90`](assemble_orbital_potential_coefficients.f90) | 读一体和各多极二体 MCP 条目，利用当前 CSF 密度矩阵将贡献收缩。按被优化轨道及积分标签分类，形成直接势、交换源和非对角项的系数列表。 |
| [`csf_density_matrix_element.f90`](csf_density_matrix_element.f90) | 用块偏移把局部 CSF 序号映射到混合数组。EOL 返回 $D_{ij}=\sum_s w_sc_{is}c_{js}$。平均层路径仅返回加权对角元。 |
| [`diagonalize_hamiltonian_block.f90`](diagonalize_hamiltonian_block.f90) | 以块对角和稀疏矩阵乘作为 Davidson 算子，首轮由小主子矩阵构造初猜，后续可复用已有本征向量。按目标根范围取出本征值和混合系数并规范相位。 |
| [`direct_angular_coefficient.f90`](direct_angular_coefficient.f90) | 零阶直接项按同壳电子对数 $q(q-1)/2$ 或两壳占据乘积计算。同壳的非零阶闭壳项按相对论 3j 系数计算。 |
| [`exchange_angular_coefficient.f90`](exchange_angular_coefficient.f90) | 当至少一个子壳层闭合时，以占据数乘积乘相对论 3j 系数平方并取负，形成对应多极阶的交换系数。 |
| [`run_scf_iterations.f90`](run_scf_iterations.f90) | 分配分块、混合和径向工作区。每轮 EOL 用当前轨道求分块哈密顿量及混合权重，更新乘子并按顺序或误差最大的轨道刷新径向函数。比较轨道误差及加权能量变化，退出后正交化并写波函数。 |
| [`solve_hamiltonian_blocks.f90`](solve_hamiltonian_blocks.f90) | 从 MCP 结构读列尾和稀疏行号，逐块代入当前径向积分，平移对角能量后求目标根。按需要跟踪和阻尼混合向量，将块平均能量及结果写入混合文件。 |
| [`update_scf_mixing.f90`](update_scf_mixing.f90) | 按等权、统计或自定义标记确定当前级权重，从密度矩阵对角收缩得到每个轨道的广义占据，并计算 $\sum_s w_s E_s$ 作为 SCF 的加权能量。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`apply_lagrange_constraints.f90`](apply_lagrange_constraints.f90) | 扫描与当前轨道有关的正交约束，按占据或代码给定的归一化权重选取乘子，将约束轨道的径向分量加入对应源项。 |
| [`apply_off_diagonal_integrals.f90`](apply_off_diagonal_integrals.f90) | 对非对角一电子积分的来源轨道计算径向导数，以能量、核势、$\kappa$ 和系数构造耦合项，更新大小分量交换源。 |
| [`build_direct_radial_potential.f90`](build_direct_radial_potential.f90) | 以核势为起点，根据直接势系数记录提取每个轨道对，逐阶计算 Hartree $Y^k$ 径向函数并按系数累加到当前轨道的 $Y(r)$ 势。 |
| [`build_exchange_radial_potential.f90`](build_exchange_radial_potential.f90) | 解包每项的径向积分轨道与源轨道，计算 Hartree 相互作用，将它与源轨道的 $P,Q$ 乘积分别累加进大小分量交换源。 |
| [`build_orbital_potential.f90`](build_orbital_potential.f90) | 依次收缩当前轨道的角系数、构造直接势及交换源，再把拉格朗日约束及非对角一电子积分并入径向方程。 |
| [`build_radial_inhomogeneous_terms.f90`](build_radial_inhomogeneous_terms.f90) | 以交换源、网格因子与步长构造 $P,Q$ 的非齐次右端并形成相邻步进需要的源数组。 |
| [`build_variation_source_term.f90`](build_variation_source_term.f90) | 以 $-Q$ 与 $P$ 分别乘半步长、Jacobian 和 $1/c$，累加相邻网格的右端项，形成随轨道能量变化的线性响应源。 |
| [`calculate_deferred_corrections.f90`](calculate_deferred_corrections.f90) | 从第三个活动网格点开始，对上一轮 $P,Q$ 施加权重为 `1,-15,40` 除以 `120` 的差分组合，得到径向步进的延迟修正，边界处置零。 |
| [`calculate_lagrange_multipliers.f90`](calculate_lagrange_multipliers.f90) | 对合格的同 $\kappa$ 轨道对，重建势和交换源后径向积分。两轨道可变时根据广义占据的相对差是否大于 `0.1`，选择势差／倒数差或势和／倒数和表达式，保存轨道对及乘子。 |
| [`check_orbital_damping.f90`](check_orbital_damping.f90) | 比较连续两次相对能量变化。启用自适应阻尼时，若乘积小于 `-1e-4` 则加强阻尼，否则将阻尼减半，并保存本次变化。 |
| [`clear_radial_inhomogeneous_terms.f90`](clear_radial_inhomogeneous_terms.f90) | 把当前大小分量的两个径向源数组置零，使下一次积分先按齐次方程推进。 |
| [`damp_orbital_update.f90`](damp_orbital_update.f90) | 按阻尼因子混合新旧大、小分量与原点参数，处理有效网格长度和符号，重新归一化并检查径向节点。 |
| [`energy_outside_bounds.f90`](energy_outside_bounds.f90) | 将候选能量与当前上下界做严格开区间比较，落在区间外或等于端点时通知更新例程采取修正。 |
| [`estimate_orbital_energy.f90`](estimate_orbital_energy.f90) | 由 $n,\kappa$ 与选取的核电荷构造氢样 Dirac 能量尺度，建立当前轨道允许的试探能量上下界。 |
| [`evaluate_radial_self_consistency.f90`](evaluate_radial_self_consistency.f90) | 在新旧轨道共有的有效网格上取 $\max_r(\lvert\Delta P\rvert+\lvert\Delta Q\rvert)$，再乘广义占据的平方根，写入当前轨道的自洽误差。 |
| [`find_least_consistent_orbital.f90`](find_least_consistent_orbital.f90) | 在可变轨道列表中查找最大的自洽误差，返回轨道索引及误差供自适应更新顺序和停止条件使用。 |
| [`improve_orbital.f90`](improve_orbital.f90) | 生成当前轨道势和差分修正，求新的径向函数与能量，按求解状态重试或作正交化。计算变化量、调整阻尼并提交更新给后续轨道。 |
| [`integrate_radial_inward.f90`](integrate_radial_inward.f90) | 组装尾部耦合方程并逐点做 Crout 消元，要求相邻尾振幅之和达到衰减容差后回代至接点。运行网格耗尽会报错。 |
| [`integrate_radial_outward.f90`](integrate_radial_outward.f90) | 从已给出的前三个原点样本开始，用变换网格的步长、$\kappa r'/r$ 及邻点源项逐步推进耦合 $P,Q$，一直到内外解接点。 |
| [`orthogonalize_same_kappa_orbitals.f90`](orthogonalize_same_kappa_orbitals.f90) | 按固定／可变轨道策略从同 $\kappa$ 轨道减去相互投影，随径向分量一起变换原点系数，再做归一化及有效网格处理。 |
| [`orthogonalize_to_more_consistent_orbitals.f90`](orthogonalize_to_more_consistent_orbitals.f90) | 选择同 $\kappa$ 且优化误差更小的参考轨道，计算重叠并减去投影。检查归一化及节点计数，才使用经修正的试探解。 |
| [`solve_radial_eigenvalue.f90`](solve_radial_eigenvalue.f90) | 用当前数值轨道分别积分交换、直接势、小分量静质量、$\kappa/r$ 和导数项，将五项相加得到轨道能量参数。 |
| [`solve_radial_orbital.f90`](solve_radial_orbital.f90) | 为试探能量建立径向系数和原点条件，先求齐次解再按需加入交换／约束源。分别由内外积分到接点，检查节点数、能量括区和接点失配，用变分源修正能量并归一化试探轨道。 |
| [`update_trial_orbital_energy.f90`](update_trial_orbital_energy.f90) | 根据当前节点数、目标节点数和能量修正移动括区。选择限制步长、区间内插或重新估计能量，并传播找不到允许能量时的状态。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`cache_mcp_files.f90`](cache_mcp_files.f90) | 先扫描一个 MCP 通道的记录数和各块边界，按实际大小分配索引／系数数组，再重新遍历填入缓存。对读入失败或损坏的记录报告 I/O 上下文。 |
| [`finish_summary.f90`](finish_summary.f90) | 输出最后一轮的能量、能级权重和组态成分，按选取的平均模式完成报告并关闭相应输出单元。 |
| [`load_mcp_records_into_memory.f90`](load_mcp_records_into_memory.f90) | 按通道编号及累计记录位置返回该积分标签和贡献数，供矩阵或势构造遍历缓存。 |
| [`open_mixing_output.f90`](open_mixing_output.f90) | 为 OL/EOL 创建混合文件并建立 `CSF--MIX` 块目录，供每轮矩阵求解写各块系数。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 将状态摘要路径连接到 SCF 文本单元，检查是否成功打开，供初始设置和结束时的能量表写入。 |
| [`print_radial_wavefunction.f90`](print_radial_wavefunction.f90) | 按列打印网格上的大、小分量及相关势的诊断样本，可调用 ASCII 绘图例程查看径向形状。 |
| [`print_scf_results.f90`](print_scf_results.f90) | 把轨道编号转换为谱学标签，列出将优化的轨道及其控制选项，供交互核对轨道选择。 |
| [`read_mcp_coefficients_from_memory.f90`](read_mcp_coefficients_from_memory.f90) | 从 `mcp.30` 缓存还原当前块的稀疏列尾偏移和行号。MPI 版本按本地列分工恢复内存布局。 |
| [`read_mcp_coefficients_into_memory.f90`](read_mcp_coefficients_into_memory.f90) | 从所选 T/V 通道缓存取一组贡献，把每条记录的 CSF 列、矩阵位置和角系数依次复制到调用方数组。 |
| [`read_mcp_structure_from_memory.f90`](read_mcp_structure_from_memory.f90) | 用对称性块在 `mcp.30` 缓存里的下标范围，将稀疏行号复制到当前块的 `matrix_partner_indices`。 |
| [`read_mcp_structure_into_memory.f90`](read_mcp_structure_into_memory.f90) | 按通道和累计贡献位置取出指定积分的 CSF 列、稀疏矩阵位置及角系数，供当前计算使用。 |
| [`write_csf_weights.f90`](write_csf_weights.f90) | 据每块 CSF、根和混合系数偏移扫过本征向量，选出幅值最大的组态分量，以块内 CSF 序号、符号和对称性输出。 |
| [`write_hamiltonian_records.f90`](write_hamiltonian_records.f90) | 遍历稀疏列的起止区间，按调试选项报告矩阵行列与元素，供核对稀疏结构和径向积分组装。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 记录核模型、离化电荷、基组规模、径向网格、优化轨道和目标能级设置，建立计算报告的输入部分。 |
| [`write_scf_energies.f90`](write_scf_energies.f90) | 对请求的根恢复块平均能量偏移，按 `2J+1` 和宇称标注级能量，打印 Hartree、cm⁻¹ 和 eV 及当前模式要求的能差。 |
| [`write_scf_radial_wavefunctions.f90`](write_scf_radial_wavefunctions.f90) | 按 `RAD--WFN` 写优化轨道和共享径向网格。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`basis_csf_parity.f90`](basis_csf_parity.f90) | 把 CSF 序号交给 RMCDHF 基组查询接口，由占据／块标签状态读取该 CSF 宇称。 |
| [`basis_csf_two_j_plus_one.f90`](basis_csf_two_j_plus_one.f90) | 从 CSF 对应的块状态取得角动量编码 $2J+1$，供统计权重和能级输出使用。 |
| [`query_available_system_memory.f90`](query_available_system_memory.f90) | 读取 `/proc/meminfo` 第二行的内存数值并换算为 GiB，用于报告缓存建立前后的系统内存读数。 |
| [`query_total_system_memory.f90`](query_total_system_memory.f90) | 读取 `/proc/meminfo` 首行的系统总内存数值，将其除以 `1048576` 后返回整数精度的 GiB，供 MCP 缓存诊断打印。 |
| [`rmcdhf_mcp_memory.f90`](rmcdhf_mcp_memory.f90) | 声明 `mcp.30` 的稀疏列界与行号，及各个 T/V 通道的标签、贡献数、列号、稀疏位置与系数数组。 |
| [`serial_mpi_context.f90`](serial_mpi_context.f90) | 声明串行内核复用的进程编号和进程数，由串行入口设为 0 和 1。 |
