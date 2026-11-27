# rwfnestimate：径向轨道初猜

`rwfnestimate` 为 MCDHF 准备初始径向轨道。它可以复用已有径向波函数，或者按
Thomas–Fermi、屏蔽氢样或自定义有效核电荷生成初猜。

程序读取 `isodata`，从 `name_header.toml` 取得轨道表和总电子数。可通过 `--header`
指定头文件路径。输出 `name.initial.wfn` 供 MCDHF 程序优化轨道。选择非默认设置时，
摘要默认写入 `name.rwfnestimate.summary`。

[MCDHF 程序](../rmcdhf/README.md) · [头文件与初猜选项](../../../docs/implemented/design-0023-rwfnestimate-header-input.md)

## 计算流程

1. 从匹配的 TOML header 装入轨道量子数、闭核数和电子数。读核输入并构造径向网格、核势与氢样屏蔽参数。[头文件输入](rwfnestimate_input.f90) · [网格和核势](read_calculation_options.f90)
2. 仅读取旧轨道时，按 $n,\kappa$ 匹配并插值到共同网格。生成初猜时可使用解析屏蔽氢样轨道，或在 Thomas–Fermi 势中求齐次 Dirac 方程。后者对试探能量分别从原点和尾部推进至接点，根据节点数及失配迭代求根。[初猜分派](generate_initial_radial_wavefunctions.f90) · [Thomas–Fermi 径向求解](solve_homogeneous_radial_equation.f90)
3. 对已初始化的轨道执行同 $\kappa$ 正交化，记录初猜来源，把波函数写成可直接加载的 `name.initial.wfn`。[正交化与写出顺序](rwfnestimate_main.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 读取轨道估算和积分的诊断位，按需打开调试输出，供径向势、网格和轨道检查使用。 |
| [`read_calculation_options.f90`](read_calculation_options.f90) | 读取 `isodata` 与可选网格参数，设置计算精度和光速，依次初始化径向网格、差分／求积权重和核 Coulomb 势。 |
| [`rwfnestimate_input.f90`](rwfnestimate_input.f90) | 解析 `--header` 或状态对应的默认路径，通过 header-only 接口读取轨道表、闭核数和电子数。由 `n,κ` 生成轨道量子数及标签，安装到共享轨道数组。 |
| [`rwfnestimate_main.f90`](rwfnestimate_main.f90) | 解析状态和头文件，读取轨道元数据及总电子数并初始化屏蔽常数。读取核和网格设置，按轨道生成初猜，最后同 κ 正交化并写径向记录。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`build_thomas_fermi_potential.f90`](build_thomas_fermi_potential.f90) | 置剩余电荷 $Z_r=Z-(N_e-1)$，以 $\sqrt[3]{N_e-1}/0.8853$ 缩放半径。用有理函数 $s(r)$ 插值屏蔽振幅，计算 $-rV_{TF}=Z_r+[-rV_{nuc}-Z_r]s^2$。 |
| [`estimate_from_saved_orbitals.f90`](estimate_from_saved_orbitals.f90) | 逐记录读取已保存的 GRASP 径向波函数，按 $(n,\kappa)$ 匹配待初始化轨道。在目标网格插值大、小分量并保存能量参数和来源标记。 |
| [`estimate_hydrogenic_orbitals.f90`](estimate_hydrogenic_orbitals.f90) | 用核电荷减去屏蔽或用菜单指定的有效电荷构造 Coulomb 场，调用解析 Dirac–Coulomb 轨道例程，保存径向分量、能量和有效范围。 |
| [`estimate_thomas_fermi_orbitals.f90`](estimate_thomas_fermi_orbitals.f90) | 逐一挑出请求 Thomas–Fermi 初猜的轨道，调用齐次径向求解器，保存成功得到的径向分量及来源。 |
| [`generate_initial_radial_wavefunctions.f90`](generate_initial_radial_wavefunctions.f90) | 初始化轨道占用标记与来源，建立 Thomas–Fermi 势，通过轨道选择菜单分批应用复用、Thomas–Fermi 或氢样初猜。每次标记已完成的轨道，并写初猜来源摘要。 |
| [`integrate_homogeneous_radial_segment.f90`](integrate_homogeneous_radial_segment.f90) | 按输入的起止下标和方向推进齐次 Dirac 大、小分量。用相邻点系数构造差分递推，重缩放中间解以控制振幅，返回积分段的径向分量。 |
| [`integrate_radial_tail.f90`](integrate_radial_tail.f90) | 由试探能量确定衰减常数及远端起点。以指数乘逆幂级数计算尾部 $P,Q$，待两分量的相对级数项均低于容差停止，再把未使用的尾点清零。 |
| [`set_hydrogenic_screening.f90`](set_hydrogenic_screening.f90) | 对闭核按顺序累加简并度 `2j+1`。当前闭壳层的屏蔽常数为此前电子数加当前简并度的一半。价层和关联轨道使用闭核电子总数。 |
| [`solve_homogeneous_radial_equation.f90`](solve_homogeneous_radial_equation.f90) | 为给定 $n,\kappa$ 建立试探能量区间，求原点和远端解，从两端积分到接点。按节点数选择括区方向，按接点失配修正能量，收敛后对 $P,Q$ 归一化并返回状态。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`open_summary_file.f90`](open_summary_file.f90) | 打开格式化的初猜摘要文件，检查 I/O 状态并关联到后续写出例程。 |
| [`print_uninitialized_orbitals.f90`](print_uninitialized_orbitals.f90) | 扫描初始化标记，把仍缺少径向函数的轨道按谱学标签列出，供交互选择下一批轨道。 |
| [`write_initial_radial_wavefunctions.f90`](write_initial_radial_wavefunctions.f90) | 按 `RAD--WFN` 写初猜轨道和共享径向网格。 |
| [`write_initial_wavefunction_summary.f90`](write_initial_wavefunction_summary.f90) | 按轨道输出初猜方法及来源说明，区分复用记录和程序生成的估计。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 写入时间、核电荷、电子／轨道计数和当前径向配置，为后续初猜来源表提供输入上下文。 |
