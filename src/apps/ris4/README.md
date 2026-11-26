## RIS4: A program for relativistic isotope shift calculations - official development version
Note: In modern Fortran and compatible with GRASP2018. First version provided by Jörgen Ekman (Malmö University) on June 7th, 2020.

### Abstract
Spectral lines from different isotopes display a small separation in energy, commonly referred to as the line isotope shift. The program ris 4 (Relativistic Isotope Shift) calculates normal and specific mass shift parameters as well as field shift electronic factors from relativistic multiconfiguration Dirac–Hartree–Fock wave functions. These quantities, together with available nuclear data, determine isotope-dependent energy shifts. Using a reformulation of the field shift, it is possible to study, in a model-independent way, the atomic energy shifts arising from changes in nuclear charge distributions, e.g. deformations.

Original publication: [Ekman et al. Comput. Phys. Comm. 235 (2019) 433](https://doi.org/10.1016/j.cpc.2018.08.017)

### RIS4 executables

`ris4` is the serial program and `ris4_mpi` distributes the matrix-element
calculation across MPI ranks. Both read `isodata`, `name_desc.parquet`,
`name_header.toml`, `name.wfn`, and either `name.mix` or `name.cimix`. They
write `name.i` for MCDHF or `name.ci` for CI. These per-state reports are the
inputs to [FICAL](#fical).

### FICAL

`fical` reads the two states' RIS4 `.i` or `.ci` reports and an `isodata` file in
the working directory. It prompts for the state names and writes
`<state1>.<state2>.fi` or `<state1>.<state2>.cfi`. It uses the unit on each
converted NMS and SMS row and reports differences in the requested GHz, MHz,
or meV units. A missing or unsupported mass-shift unit is an input error.
The linear-dependence discrepancy is reported as `NaN` when its normalization denominator is zero.

The Suite builds and installs `fical` as a separate executable linked against
`core`. After configuring Suite, build the target from the repository root.
Replace `build-gf16` with your configured build directory. The executable is
written to that directory's `bin/`. [Suite build setup](../../../README.md#build)

```fish
cmake --build build-gf16 --target fical
```

### 计算流程

1. `ris4` 和 `ris4_mpi` 在给定基组、混合向量及径向轨道上预计算轨道对的接触密度、径向矩、NMS 和反冲径向积分。以一体、二体角系数与径向表相乘，再与态的混合系数收缩，得到各级电子因子。MPI 版本对 CSF 对分工并归约各级结果。[串行位移计算](calculate_isotope_shifts.f90) · [MPI 位移计算](calculate_isotope_shifts_mpi.f90)
2. Seltzer 路径额外累计近核密度随 $r$ 的剖面，用 $\rho(r)=p_0+p_2r^2+p_4r^4+p_6r^6$ 拟合，并把各项换算成与核径向矩差对应的电子因子。报告同时记录质量位移及场位移所需的数据。[Seltzer 密度收缩](calculate_seltzer_field_shift_factors.f90) · [电子密度拟合](fit_electron_density.f90)
3. `fical` 读两态报告，统一每行的质量位移单位，按能级配对取得上下态因子差。质量项结合核质量倒数差，场位移项按选用的径向矩或半径模型计算，最终报告分项及线位移。[谱线同位素位移](isotope_shift_transition_main.f90)

### 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

#### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 初始化诊断位，读取允许的非默认诊断选择，按需连接 `.dbg` 单元，供输入、积分与结果检查使用。 |
| [`prepare_isotope_calculation.f90`](prepare_isotope_calculation.f90) | 解析核电荷、核质量、核分布及非默认轨道设置，生成网格、差分表和核势，再读取径向波函数供所有位移算符使用。 |
| [`read_mixing_blocks.f90`](read_mixing_blocks.f90) | 读取状态对应的 `CSF--MIX`，逐块提取根编号、$J$、宇称、绝对能量和本征向量，将其按 CSF 偏移铺到共享数组，供各算符按态收缩。 |
| [`ris4_main.f90`](ris4_main.f90) | 选择 MCDHF／CI，装入状态基组、核和径向数据，再读各级混合向量，调用 RIS4 的电子因子计算并输出每态报告。 |
| [`ris4_mpi_main.f90`](ris4_mpi_main.f90) | 由 rank 0 解析设置并打开报告，集体读入描述符基组。载入径向、核和混合数据后调用 MPI 位移计算。 |

#### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`first_spin_angular_factor.f90`](first_spin_angular_factor.f90) | 用 bra 的反号 $\kappa$ 对应轨道角动量及双方的 $j$ 构造 Racah 重耦合，返回第一个 SMS 自旋约化因子。 |
| [`second_spin_angular_factor.f90`](second_spin_angular_factor.f90) | 用 ket 的反号 $\kappa$ 对应轨道角动量构造另一重耦合，返回第二个 SMS 自旋约化因子。 |

#### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`field_shift_radial_integral.f90`](field_shift_radial_integral.f90) | 在径向网格下标 10、12、14 计算轨道对密度 $(P_aP_b+Q_aQ_b)/(4\pi r^2)$。用偶次插值外推到 $r=0$。 |
| [`field_shift_radial_integral_vector.f90`](field_shift_radial_integral_vector.f90) | 把同一个轨道对的 $(P_aP_b+Q_aQ_b)/(4\pi r^2)$ 逐点写入近核网格数组，从第二个网格点开始，供 Seltzer 拟合路径收缩。 |
| [`interpolate_electron_density.f90`](interpolate_electron_density.f90) | 解三点的偶次多项式插值问题，用常数项估计原点电子密度。 |
| [`normal_mass_shift_radial_integral.f90`](normal_mass_shift_radial_integral.f90) | 差分求 $P,Q$ 的径向导数，积分导数乘积和由 $l,\kappa$ 决定的径向项，再按相对论核电场项增加修正，返回领先项与完整 NMS 积分。 |
| [`specific_mass_shift_radial_integral_2.f90`](specific_mass_shift_radial_integral_2.f90) | 调用两个自旋角因子，按其权重组合 $P_aQ_b,Q_aP_b$ 与核反冲径向因子，再求积得到第二项。 |
| [`specific_mass_shift_radial_integral_3.f90`](specific_mass_shift_radial_integral_3.f90) | 构造小分量乘积与径向核因子的被积函数，在共同有效区间积分得到第三项。 |
| [`vinti_radial_integral.f90`](vinti_radial_integral.f90) | 对径向分量求导，用两轨道 $\kappa$ 推导的角因子构造 ${1}/{r}$ 项，在共同有效网格上积分出 Vinti 因子。 |

#### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`calculate_field_shift_factors.f90`](calculate_field_shift_factors.f90) | 用一电子标量角系数与轨道接触密度、普通质量位移和径向矩相乘，生成 CSF 对的矩阵元，再乘各态混合系数并累加至能级。 |
| [`calculate_field_shift_factors_mpi.f90`](calculate_field_shift_factors_mpi.f90) | 按 rank 分配 CSF 对并预筛子壳层占据差，复用每对 CSF 的混合乘积来累加接触密度、NMS 和径向矩，按选项保存本地角系数。 |
| [`calculate_isotope_shifts.f90`](calculate_isotope_shifts.f90) | 预制各轨道对的接触密度、径向矩、普通质量位移和反冲积分表。按选项重算或读取一／二体角系数，收缩为各级 NMS、SMS 和场位移因子，必要时拟合近核密度并写单位换算结果。 |
| [`calculate_isotope_shifts_mpi.f90`](calculate_isotope_shifts_mpi.f90) | 协调并广播轨道积分与角文件模式，区分串行缓存和 rank 记录。按进程分配一体、二体矩阵元，归约位移和电子密度，再由 rank 0 计算拟合及报告。 |
| [`calculate_seltzer_field_shift_factors.f90`](calculate_seltzer_field_shift_factors.f90) | 在标量一体矩阵元之外，把轨道密度剖面按角系数累加成 CSF 对的近核径向密度，再按各级混合系数收缩到 `level_electron_density_profiles`。 |
| [`calculate_seltzer_field_shift_factors_mpi.f90`](calculate_seltzer_field_shift_factors_mpi.f90) | 在分配的 CSF 对上累加各级近核网格密度、接触密度、NMS 和径向矩，形成供统一拟合的 rank 局部剖面。 |
| [`calculate_specific_mass_shift.f90`](calculate_specific_mass_shift.f90) | 对 CSF 对计算二体角系数，按打包轨道四元组提取 Vinti 及相对论修正积分。组合出首项和总 SMS 矩阵元后与各态混合系数收缩。 |
| [`calculate_specific_mass_shift_mpi.f90`](calculate_specific_mass_shift_mpi.f90) | 预筛占据差并循环派发 CSF 列，对本 rank 的二体项应用角系数和反冲积分，积累各态的首项及完整 SMS 贡献。 |
| [`fit_electron_density.f90`](fit_electron_density.f90) | 选取可信的近核网格，先用三点外推确定 $\rho(0)$，再用加权最小二乘拟合剩余 $r^2,r^4,r^6$ 项。将密度系数换算成 $F_0,F_2,F_4,F_6$，并计算拟合偏差与有效高阶贡献。 |
| [`isotope_shift_transition_main.f90`](isotope_shift_transition_main.f90) | `fical` 读两份态报告和同位素数据，逐行识别 NMS／SMS 单位并统一换算，按能量及能级配对取上下态电子因子差。由核质量倒数差求质量位移，以核径向矩差或输入半径模型求场位移。输出合计和线性偏离指标。 |

#### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`open_isotope_shift_angular_cache.f90`](open_isotope_shift_angular_cache.f90) | 按算符类别找 `.IOB` 或 `.ITB`，先探测底层流可读性，再检查最后一条顺序记录是否为 `-1`。完整缓存回卷到开头。缺失或未写完的缓存可重新计算。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 按状态和混合文件类型为 `.i` 或 `.ci` 建立报告输出单元，交给电子因子例程写入。 |
| [`read_field_shift_factors.f90`](read_field_shift_factors.f90) | 读取缓存的 CSF 对和单电子角系数。重建接触密度、NMS 和径向矩的矩阵元，再用当前混合系数得到能级值。 |
| [`read_field_shift_factors_mpi.f90`](read_field_shift_factors_mpi.f90) | 按记录所属 rank 取缓存的一体角系数，以对应轨道积分组装局部的 $\rho(0)$、NMS 和径向矩贡献。 |
| [`read_seltzer_field_shift_factors.f90`](read_seltzer_field_shift_factors.f90) | 复用缓存一体角系数并同时重建每个 CSF 对的径向密度剖面，收缩到能级后交给拟合。 |
| [`read_seltzer_field_shift_factors_mpi.f90`](read_seltzer_field_shift_factors_mpi.f90) | 用 rank 记录筛选缓存条目，重建本地一体矩阵元和近核密度剖面，按态累加等待归约。 |
| [`read_specific_mass_shift.f90`](read_specific_mass_shift.f90) | 解包缓存的二体角标签得到轨道四元组，代入 Vinti 表与修正表，再按各态向量重算首项和总 SMS。 |
| [`read_specific_mass_shift_mpi.f90`](read_specific_mass_shift_mpi.f90) | 按 rank 所有权筛选二体角缓存，解包轨道标签，用共同径向表重算各能级的局部 SMS 贡献。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 向 RIS4 报告写电子、轨道和 CSF 规模等输入摘要。 |
| [`write_level_energies.f90`](write_level_energies.f90) | 按所选能级写能量和可选相邻能差或相对最低态的能差，分别用 Hartree、cm⁻¹ 与 eV 表示。 |

#### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`isotope_shift_options.f90`](isotope_shift_options.f90) | 声明各能级的接触密度、径向矩，以及首项和完整 NMS／SMS 的结果数组，供角向收缩及报告共用。 |
