# rtransition_phase：带符号跃迁矩阵元

`rtransition_phase` 计算电、磁多极跃迁，并在主报告中额外写出带符号的跃迁矩阵元。
其大小为线强度的平方根，符号取自计算得到的跃迁振幅。电多极跃迁分别输出
Coulomb 和 Babushkin 规范的值。

输入与 `rtransition` 相同，包括 `isodata`、两态的描述符文件对及双正交轨道和混合系数。
带符号的值写入 `initial.final.mcdhf.transition.phase` 或 `initial.final.ci.transition.phase`。
两态均有配套 LSJ 标签时另写 `.transition.phase.lsj`。

[基本跃迁程序与输入文件](../rtransition/README.md) · [矩阵元输出实现](write_transition_probability.f90)

双态跃迁读取与 `rtransition` 共用的有序态对目录。同态跃迁直接读取对应态的原始轨道和混合流。 [路径与校验规则](../../../docs/reviews/evidence-0033-output-file-formats.md#有序态对交接与跃迁报告)

## 计算流程

1. 读取初末态描述符与双正交波函数，装配多极角记录和轨道对径向积分。对选定的能级对分别收缩得到电跃迁两规范或磁跃迁的带符号振幅。[入口](rtransition_phase_main.f90) · [振幅计算](calculate_oscillator_strengths.f90)
2. 对每个写出的结果，先由振幅平方得到非负线强 $S$，然后用该规范的原振幅符号构造 $\operatorname{sgn}(M)\sqrt{S}$。主报告在跃迁率和线强外记录这项矩阵元。[带相位的输出实现](write_transition_probability.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`load_final_lsj_labels.f90`](load_final_lsj_labels.f90) | 将末态 LSJ 标签文件的标签和级记录装入末态标签数组，供配对时使用。 |
| [`load_final_transition_basis.f90`](load_final_transition_basis.f90) | 加载末态描述符，保存电子数、闭核、块前缀和轨道编码。安装用于合并的另一套 CSF 数据。 |
| [`load_initial_lsj_labels.f90`](load_initial_lsj_labels.f90) | 打开初态 LSJ 标签文件并读入根编号、LSJ 标签和相关能量字段，供可选谱线输出匹配。 |
| [`load_initial_transition_basis.f90`](load_initial_transition_basis.f90) | 加载初态描述符，保存其闭核、块前缀、轨道、占据和累计耦合数组，以便载入另一态后仍能访问原始 CSF。 |
| [`map_merged_orbitals.f90`](map_merged_orbitals.f90) | 对合并后的轨道表逐项匹配原初、末态轨道，建立公共索引到两套径向数组槽位的映射。 |
| [`merge_transition_basis_records.f90`](merge_transition_basis_records.f90) | 调用共享基组合并器，在内存中统一轨道、重映射占据与耦合。无法合并时打印两态名和错误原因并终止。 |
| [`merge_transition_csf_bases.f90`](merge_transition_csf_bases.f90) | 依次加载两套描述符基组并连接基组归并例程，让后续张量角计算按共同轨道顺序访问两套态。 |
| [`read_final_radial_wavefunctions.f90`](read_final_radial_wavefunctions.f90) | 读取末态 轨道文件 的径向数组及量子数，以末态索引保存匹配记录供跨态积分使用。 |
| [`read_initial_radial_wavefunctions.f90`](read_initial_radial_wavefunctions.f90) | 打开初态 轨道文件，验证二进制标记并按 $n,\kappa$ 匹配网格记录，载入初态轨道的分量、能量和原点参数。 |
| [`read_multipole_options.f90`](read_multipole_options.f90) | 把 `E1`、`M2` 等词元解析成同时编码阶数和宇称的整数，拒绝不可识别的输入，并读取 M1/E2 是否限于异 $J$ 态。 |
| [`calculation_state_files.f90`](../../lib/core/calculation_state_files.f90) | 读取初态与末态状态名，按顺序交给文件加载和跨态输出路径解析。 |
| [`read_transition_mixing_coefficients.f90`](read_transition_mixing_coefficients.f90) | 读取已定位的初、末态混合块，分别提取各级编号、对称性、平均能量、本征值和 CSF 向量。 |
| [`read_transition_options.f90`](read_transition_options.f90) | 设置核和径向网格、选择电／磁多极、能量单位、阈值和输出模式，装入两态双正交轨道并检查跨态径向重叠。 |
| [`rtransition_phase_main.f90`](rtransition_phase_main.f90) | 沿串行入口的初末态合并、径向积分及混合系数收缩路径，使用带符号矩阵元的输出例程生成相位版主报告。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`sort_transition_angular_coefficients.f90`](sort_transition_angular_coefficients.f90) | 生成排序索引和各积分标签的起止位置，对相同初末轨道对归组，保存每组的 CSF 对与角因子记录。 |
| [`write_transition_angular_coefficients.f90`](write_transition_angular_coefficients.f90) | 遍历指定初末对称块的 CSF 对，按所选多极的秩、宇称生成单电子张量角系数。按打包轨道对排序保存供积分收缩。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`build_transition_bessel_functions.f90`](build_transition_bessel_functions.f90) | 根据跃迁频率、多极阶和径向长度准备球 Bessel 函数表，对不同自变量范围用适合的级数或递推，保存供规范积分使用的阶次。 |
| [`calculate_biorthogonal_overlap.f90`](calculate_biorthogonal_overlap.f90) | 在初、末轨道共用的径向区间求大、小分量乘积之和的积分，给出双正交准备完成后的检查矩阵。 |
| [`single_particle_transition_element.f90`](single_particle_transition_element.f90) | 由合并轨道索引定位初末 $P,Q$，与需要的球 Bessel 阶次及相对论角因子组合后求积。电多极返回 Coulomb 和 Babushkin 振幅，磁多极返回相应磁振幅。 |

### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`calculate_oscillator_strengths.f90`](calculate_oscillator_strengths.f90) | 逐多极和初末对称块读取混合与角记录，应用能级选择规则，恢复块偏移求跃迁能量。用该频率重建径向核，收缩两态向量得到振幅，再输出跃迁率和线强度。 |
| [`write_transition_csf_weights.f90`](write_transition_csf_weights.f90) | 按轨道积分分组求单电子径向矩阵元，把每条角系数乘以初、末态相应 CSF 混合系数的乘积，累加到级间振幅，按需打印主要贡献。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`check_transition_mixing_coefficients.f90`](check_transition_mixing_coefficients.f90) | 将已读进的两态维数、恢复偏移后的本征值、级和 CSF 系数完整打印供诊断。 |
| [`open_transition_angular_cache.f90`](open_transition_angular_cache.f90) | 核对块头和末尾完成记录中的初末态描述符摘要、多极代码、J 筛选及进程布局。旧格式、截断或身份不符的角缓存会重算。[共享格式](../../lib/core/transition_angular_cache.f90) |
| [`read_transition_angular_coefficients.f90`](read_transition_angular_coefficients.f90) | 读取指定块／多极的积分数、打包轨道对、每组的贡献边界及系数，检查并扩充工作数组后交给振幅计算。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 由初、末态名及混合来源构造并打开跃迁报告，启用标签时建立对应的 LSJ 报告。 |
| [`write_lsj_transition_probability.f90`](write_lsj_transition_probability.f90) | 把级编号映射为初末 LSJ 标签，配合频率或波长及两规范的强度，写出可选的带标签谱线记录。 |
| [`write_transition_level_energies.f90`](write_transition_level_energies.f90) | 把本征值加回块平均位移后列出两态各级总能量，按 Hartree、cm⁻¹ 与 eV 生成可比较的能级数据。 |
| [`write_transition_probability.f90`](write_transition_probability.f90) | 用同一振幅和频率生成 `A`、`gf`、`S`，主报告的每个适用行还写 $\operatorname{sgn}(M)\sqrt{S}$。电跃迁分别保留两个规范的符号，磁跃迁保留磁振幅的符号。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`resize_transition_coefficient_arrays.f90`](resize_transition_coefficient_arrays.f90) | 按所选操作为角系数和相关 CSF 记录数组分配、扩容或释放，维持同一容量边界。 |
| [`resize_transition_integral_arrays.f90`](resize_transition_integral_arrays.f90) | 为不同积分的轨道标签、起止偏移等数组同步扩容或释放，保留已有分组记录供继续读取。 |
| [`resize_transition_label_arrays.f90`](resize_transition_label_arrays.f90) | 管理角系数排序时的打包标签及关联索引缓冲容量。 |
| [`resize_transition_matrix_arrays.f90`](resize_transition_matrix_arrays.f90) | 同步分配或调整跃迁矩阵贡献数组的空间，供读取角记录和收缩振幅使用。 |
| [`right_csf_intermediate_coupling.f90`](right_csf_intermediate_coupling.f90) | 通过当前跨态配对状态的中间耦合查询，按子壳层位置与 CSF 序号返回耦合编码。 |
| [`right_csf_parity.f90`](right_csf_parity.f90) | 按 CSF 序号返回对应工作数组中的宇称标记，供跃迁角动量选择。 |
| [`right_csf_subshell_occupation.f90`](right_csf_subshell_occupation.f90) | 以合并轨道和 CSF 序号读取指定子壳层占据，为一电子算符判断可连接的 CSF 对。 |
| [`right_csf_subshell_quantum_number.f90`](right_csf_subshell_quantum_number.f90) | 按子壳层、量子数槽及 CSF 序号从工作基组读取相应的占据态量子数。 |
| [`right_csf_two_j_plus_one.f90`](right_csf_two_j_plus_one.f90) | 获取工作基组中对应 CSF 的 $2J+1$ 编码，为块与多极三角选择提供输入。 |
