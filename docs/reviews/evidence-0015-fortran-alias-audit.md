# Suite 的 Fortran 别名与指针关联审计

2026-11-20，源码基准为 `e8a15a9`，统计包含命名清理前的工作树修改。前两批清理消除了导入重命名，现行接口见末节的旧名迁移表。

Fortran 源码中的 **1,014 处 `=>` 有 959 处是 `USE, ONLY` 导入别名，占 94.6%**。增长主要来自共享库命名迁移，重复映射最集中在自旋–角向库。具体发现包括 RCI 串行和 MPI 共有的导入冲突、角向过程的无用导入，以及 Davidson 的 7 处无用指针关联。

## 统计口径与版本变化

扫描 `src/` 和 `test/` 中 1,062 个现存的 `.f90` 和 `.F90` 文件。统计按符号出现次数计，合并自由格式续行，排除注释和字符串。含预处理分支的文件按源码计。原始文本匹配 1,030 次，其中 16 次在注释或字符串中。Rust 的 `match` 分支和 `.f90_PROFILING` 备用文件在 Fortran 统计范围之外。

| Fortran 用法 | 生产源码 | 测试 | 合计 |
| --- | ---: | ---: | ---: |
| `use m, only: local => original` | 844 | 115 | 959 |
| `ptr => target` 数据或过程指针关联 | 26 | 7 | 33 |
| `pointer :: p => null()` 空关联初始化 | 20 | 0 | 20 |
| `associate (entry => entries(i))` 块内名称关联 | 2 | 0 | 2 |
| 合计 | 892 | 122 | 1,014 |

| 版本 | 导入别名 | 全部语法用例 |
| --- | ---: | ---: |
| `58fc01a`，检查时本地 `origin/main` 指向的版本 | 414 | 469 |
| `2d4c721`，共享库命名迁移前 | 414 | 469 |
| `7097193`，共享库命名与接口迁移 | 906 | 961 |
| `e8a15a9`，检查时的 `HEAD` | 961 | 1,016 |
| 检查时的工作树 | 959 | 1,014 |

早期的 `4c17db6` 先在定义端改名。例如，`nw` 成为 `orbital_count` 后，部分测试通过 `nw => orbital_count` 继续使用原名。该提交同时移除了 CSFG。在前后均存在的 1,901 个 Fortran 文件中，含 `ONLY` 和不含 `ONLY` 的导入别名合计由 201 处增至 296 处。随后 `7097193` 在角向过程的使用端引入描述性名称，全仓净增 492 处导入别名。表中各版本的指针关联、空关联初始化和 `ASSOCIATE` 次数一致。

## 角向库的重复映射与共享状态

`src/lib/spin_angular/` 有 511 处导入别名，占导入总数的 53.3%。`src/lib/core/` 有 115 处，RCI 串行和 MPI 各有 39 处。按被导入的模块统计，前三个模块合计占 51.5%。

| 被导入的模块 | 导入别名 | 不同的映射 | 示例 |
| --- | ---: | ---: | --- |
| `jj_tensor_coupling_tables` | 240 | 16 | `ket_shell_values_1 => bd1` |
| `csf_pair_angular_state` | 175 | 10 | `bra_occupations => nq1` |
| `angular_symmetry_labels` | 79 | 3 | `right_justified_j_labels => jlbr` |

`USE, ONLY` 的左侧是当前作用域访问右侧模块实体时使用的名字。例如，`setqna` 通过 `bra_occupations => nq1` 写入当前 CSF 对的 bra 占据数，`onescalar2` 读取同一数组并计算转移相位。`perko2` 同样通过模块数组把活动子壳层的数据传给角系数例程。这些名称明确了物理角色，共享数组的读写顺序仍由调用链决定。[CSF 对工作区](../../src/lib/grasp_shared/csf_pair_angular_state.f90) · [CSF 对量子数准备](../../src/lib/core/csf_pair_quantum_numbers.f90) · [标量转移系数](../../src/lib/spin_angular/one_electron_scalar_between_subshells.f90) · [活动子壳层准备](../../src/lib/spin_angular/jj_active_subshell_states.f90)

`jj_tensor_coupling_tables` 的 16 个映射在 15 个过程中整组导入。按过程检查 `USE` 列表之外的标识符引用，其中 95 处导入无引用。`onescalar2` 单个过程的 27 处导入别名中有 16 处无引用，包括第三、第四个活动子壳层的数据和编码。该过程调用 `perko2` 时指定两个活动子壳层，随后也只使用前两个子壳层的张量数据。[张量数组定义](../../src/lib/grasp_shared/jj_tensor_coupling_tables.f90) · [onescalar2](../../src/lib/spin_angular/one_electron_scalar_between_subshells.f90)

## RCI 串行和 MPI 的同名导入冲突

两个版本的 `prepare_correction_integral_tables` 都用下列本地名称导入两张不同的积分表。本地名称中带有 `vacuum_polarization`，其中一组实际来自库仑积分表。

| 本地名称 | `coulomb_integral_table` | `vacuum_polarization_integral_table` |
| --- | --- | --- |
| `vacuum_polarization_integral_indices` | `indoei` | `indvpi` |
| `vacuum_polarization_integral_values` | `valoei` | `valvpi` |

清理前两版过程只在导入列表中使用这两个名称。使用当时的模块和 GNU Fortran 15.3.0 验证时，仅导入的例程通过了语法检查，引用两个名称的例程分别报出 `ambiguous reference`。[串行 RCI](../../src/apps/rci/prepare_correction_integral_tables.f90) · [MPI RCI](../../src/apps/rci_mpi/prepare_correction_integral_tables.f90)

在工作区根目录下创建 `temp/tests/fortran-arrow-rci/`，将以下程序保存为其中的 `rci_import_probe.f90`。

```fortran
program rci_import_probe
    use coulomb_integral_table, only: vacuum_polarization_integral_indices => indoei, &
        vacuum_polarization_integral_values => valoei
    use vacuum_polarization_integral_table, only: vacuum_polarization_integral_indices => indvpi, &
        vacuum_polarization_integral_values => valvpi
    implicit none (type, external)
#ifdef REFERENCE_ALIAS
    print *, size(vacuum_polarization_integral_indices)
    print *, size(vacuum_polarization_integral_values)
#endif
end program rci_import_probe
```

从工作区根目录提取清理前的两个模块。第一条编译命令只检查导入，第二条启用两处引用，预期报两处歧义错误。

```fish
set -l audit_dir temp/tests/fortran-arrow-rci
set -l shared_src src/lib/grasp_shared
for module_name in coulomb_integral_table vacuum_polarization_integral_table
    git -C grasp-atomic-suite show "e8a15a9:$shared_src/$module_name.f90" > "$audit_dir/$module_name.f90"
end
/usr/bin/gfortran-15 -std=f2018 -cpp -fsyntax-only -J "$audit_dir" -I "$audit_dir" \
    "$audit_dir/coulomb_integral_table.f90" "$audit_dir/vacuum_polarization_integral_table.f90" \
    "$audit_dir/rci_import_probe.f90"
/usr/bin/gfortran-15 -std=f2018 -cpp -DREFERENCE_ALIAS -fsyntax-only -J "$audit_dir" -I "$audit_dir" \
    "$audit_dir/rci_import_probe.f90"
```

## 指针与块内名称关联

- RMCDHF 先分配新占据数组并读取全部描述符块，成功后才释放旧数组，通过 `csf_subshell_occupations => new_occupations` 让模块指针指向新数组。共享的 `ralloc` 例程则先把需保留的数据复制到新数组，再释放旧存储并切换指针。两类关联都承担数组替换。[RMCDHF 描述符输入](../../src/lib/core/rmcdhf_csf_input.f90) · [数组存储管理](../../src/lib/grasp_shared/array_memory_management.f90)
- Davidson 的 `iniest` 把 7 个局部指针关联到可分配工作数组，之后都直接用可分配数组执行 LAPACK 调用或访问数据。这 7 个指针在建立关联后均无引用。[稀疏 Davidson 初始估计](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90)
- `rlevels` 的 `associate (entry => entries(index_in_table))` 在块内简写数组元素。对 `entry%configuration` 的赋值更新原数组中的能级记录。[rlevels](../../src/apps/rlevels/rlevels.f90)

## 审计时确定的清理范围

1. RCI 两个版本可移除造成冲突且无引用的导入项。5 处同名映射也可直接简写为普通 `ONLY` 项，其中 `orbital_count => orbital_count` 有 2 处，`electron_count => electron_count` 有 3 处。[MPI 角向输入](../../src/apps/rangular_mpi/rangular_mpi_input.f90) · [描述符读取](../../src/lib/core/csf_descriptor_reader.f90) · [描述符对照测试](../../test/descriptor_oracle/csf_descriptor_verify.f90)
2. 角向过程可先移除无引用的导入项，再评估把表中反复出现的物理名称移到共享模块的定义端。具备多种用途的工作数组需要按上下文判断。`radial_quadrature_workspace` 中的 `ta` 在求积过程里表示加权被积函数，在其他过程中还会暂存系数或径向导数。[径向求积工作区](../../src/lib/grasp_shared/radial_quadrature_workspace.f90) · [径向源项](../../src/lib/core/radial_origin_solution.f90)
3. Davidson 的局部清理可针对上述 7 个指针及其关联。涉及实际替换存储的指针需要核对分配、失败和释放路径。共享容量模块中的 6 个旧名称是兼容接口，包含 `nnnw => max_orbitals`。[共享容量重导出](../../src/lib/grasp_shared/grasp_parameters.f90) · [共同参数契约](../implemented/design-0003-shared-grid-parameters.md)

以上结论来自清理前的静态检查与 RCI 最小编译复现。

## 清理后的命名与验证

2026-11-20，第一批清理以修改开始时的工作树为对照，887 份生产 Fortran 源文件有改动，移除了 3 份仅含废弃状态声明的应用模块。共享模块直接声明 `bra_occupations`、`ket_shell_codes_1`、`requested_level_indices` 等名称，计算过程通过 `USE, ONLY` 引用同名实体。`spin_angular` 中原有的 511 处导入别名全部消除。过程调用者可查阅 [186 项共享库过程名称对照](../implemented/design-0007-shared-library-procedure-names.md)。

第二批清理处理了余下的 34 处导入别名。当时的扫描覆盖 `src/` 和 `test/` 中全部 1,060 份 Fortran 源文件，使用与初次统计相同的续行、注释和字符串处理口径。

| Fortran 用法 | 清理前 | 第一批后 | 最终生产源码 | 最终测试 | 最终合计 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `USE` 导入重命名 | 959 | 34 | 0 | 0 | 0 |
| 指针关联赋值 | 33 | 26 | 19 | 7 | 26 |
| 指针空关联初始化 | 20 | 20 | 20 | 0 | 20 |
| `ASSOCIATE` | 2 | 2 | 2 | 0 | 2 |

| 第二批处理的导入别名 | 数量 | 现行实施方式 |
| --- | ---: | --- |
| 径向工作数组 | 14 | 各过程直接使用共享数组，注释说明该处存放的物理量。例如 Slater 积分用工作数组 B 接收 Hartree Y，用 A 形成带权被积函数。 |
| 容量定义和容量测试 | 12 | `grasp_parameters` 直接导出六个语义容量常量。容量测试核对这些数值、派生关系和实际数组维度。 |
| 节点计数与 CSF 字节打包 | 2 | 两个公开模块直接定义 `count_radial_nodes` 和 `pack_csf_subshell_byte`。生产及测试均使用正式过程。 |
| Breit 角系数缓冲 | 2 | 外部 `collect_breit_angular_coefficient` 负责链接分派，`store_breit_angular_coefficient` 提供默认存储。RCI 两个可执行文件分别提供自己的收集实现。 |
| RDENSITY 基组加载 | 2 | 串行与 MPI 的 `density_parallel` 模块各自实现 `load_density_basis`，保留各自原有的关键字形参并调用相应加载器。 |
| 描述符测试的输出单元 | 2 | 直接使用 `iso_fortran_env` 的 `output_unit`。RMCDHF 夹具给占据数二进制文件单独使用 `occupations_unit`。 |

正式过程支持关键字、过程实参和过程指针。Breit 测试分别验证新的外部收集入口和直接存储路径。[公开过程测试](../../test/f2018/library_modules.f90) · [Breit 链接替换测试](../../test/f2018/breit_procedure_hooks.f90)

RCI 两版均已移除重名的积分表导入。Davidson 直接使用可分配工作数组。RMCDHF 在全部描述符块读入成功后切换占据数组的指针关联。[径向求积工作区](../../src/lib/grasp_shared/radial_quadrature_workspace.f90) · [RCI 积分表准备](../../src/apps/rci/prepare_correction_integral_tables.f90) · [Davidson 初始估计](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90) · [RMCDHF 描述符读取](../../src/lib/core/rmcdhf_csf_input.f90)

应用命名按公式和数据流核对。例如，RIS4 的旧 `dint2` 在生产过程中接收含相对论修正的总正常质量位移积分，旧 `dint7` 接收首项，现名分别为 `total_nms_integrals` 和 `leading_nms_integrals`。两者的区别可以从径向例程的赋值以及输出报告中“总量减首项”的修正量复核。RMCDHF 中原先分时存放权重和、占据数和、加权能量的累加变量拆成了三个变量。原子角动量记号、标准线性代数例程、输入字段及外部链接钩子按各自约定命名。[RIS4 径向公式](../../src/apps/isotope/normal_mass_shift_radial_integral.f90) · [RIS4 积分表与输出](../../src/apps/isotope/calculate_isotope_shifts.f90) · [RMCDHF 混合系数更新](../../src/apps/rmcdhf/update_scf_mixing.f90)

| 第一批验证 | 结果 |
| --- | --- |
| GNU Fortran 15.3 Check，构建串行与 MPI 程序，启用 jj2lsj OpenMP | 116/116 项 CTest 通过 |
| GNU Fortran 15.3 Release，关闭 MPI 与 jj2lsj OpenMP | 76/76 项 CTest 通过 |
| 同一编译器和 Check 配置，清理前后逐字节比较 | 63 次文件对照一致。计算链覆盖串行及双进程 MPI 的轨道、混合系数和 MCP，性质计算覆盖双正交变换、跃迁、Landé 因子、超精细及同位素位移，RCI 包含 Coulomb 与 Breit 两种设置 |
| 过程体静态对照 | 核对 1,113 个作用域中的 49,058 条可执行语句。按名称迁移和删除的 7 个无引用指针关联对齐后，发现的一处字面量差异是 Hartree 核张量阶数越界时的诊断文本，已改用实际参数名和过程名 |
| RMCDHF 来源核查 | 212 个改动文件连接到前一份来源快照，现存 250 个受保护文件通过完整性核查及相应篡改回归。[命名快照](../reference/rmcdhf_semantic_names_snapshot.json) |

第二批完成后，GNU Fortran 15.3 的 Check 配置通过 117/117 项 CTest，纯串行 Release 通过 77/77 项。网格补丁器通过 37 项测试，包括历史导入布局以及 590、1990 两种容量的编译验证。使用安装目录中的模块与静态库编译独立调用程序，也通过了具体过程、关键字和 Breit 外部链接替换检查。

以仍有 34 处导入别名的 Check 二进制为对照，串行及双进程 MPI 的 SCF 流程分别得到 16 份逐字节一致的文件。RIS4 两种设置共比对 6 份文件，Coulomb 与 Breit 两种 RCI 设置共比对 2 份，串行与 MPI 性质计算分别比对 12 和 11 份，合计 63 次文件对照全部一致。

外部 Fortran 调用方需要用相应版本的 `.mod` 和静态库重新编译，并按对照表更新已改名的导入。上述验证覆盖接口、执行路径和指定夹具的输出。变量名的物理含义以对应公式、数组定义及调用方为依据。

## 旧命名入口的迁移

第三批移除了先前为旧调用方保留的常量、过程、关键字和模块。下表左列用于查找需要迁移的源代码。

| 旧接口 | 当前调用方式 |
| --- | --- |
| `KEYORB`、`NNNP`、`NNN1` | 分别使用 `orbital_index_radix`、`max_radial_points`、`radial_workspace_capacity` |
| `NNNW`、`NNNWM1`、`NNNWM2` | 分别使用 `max_orbitals`、`max_orbitals_minus_one`、`max_orbitals_minus_two` |
| `count`、`pack` | 从 `radial_node_count` 导入 `count_radial_nodes`，从 `csf_byte_packing` 导入 `pack_csf_subshell_byte`。字节打包的关键字为 `unpacked_value`、`byte_index`、`packed_values` |
| `gdvd`、`rkco_gg` | 使用 `solve_davidson_eigenpairs` 和 `calculate_two_electron_angular_coefficients`。Davidson 调用删除原位置实参 `orthogonality_tolerance`，保留能量、系数和残差三个实际收敛判据 |
| `talk`、`buffer_breit_angular_coefficient`、`talk_fallback` | 需要链接时替换的调用使用 `collect_breit_angular_coefficient`，直接写入默认缓冲区使用 `store_breit_angular_coefficient`。外部替换实现也须改名并重新链接 |
| RCI 的 `dspevx_interface` | 直接从 `blas_lapack_interfaces` 导入 LAPACK 的 `dspevx` |
| Biotransform 矩阵求逆的旧关键字 | 串行与 MPI 入口使用 `matrix`、`inverse_workspace`、`matrix_order`、`determinant`、`singularity_tolerance`、`inversion_status`、`leading_dimension` |

`mpi_state` 中仅供早期接口使用的 `ierrtotal` 已移除。网格工具针对 `atomic-suite` 校验六个规范容量名。工具处理独立 GRASP2018 源码时继续使用该项目自身的命名。两种布局共用的 CLI 和 TOML 参数 `nnnp` 在 Suite 中修改 `max_radial_points`。[容量定义与网格配置](../implemented/design-0003-shared-grid-parameters.md)

另有 71 个过程曾将中央轨道打包进制定义为同值局部常量，其中 70 处使用 `orbital_label_radix`，1 处使用 `coefficient_label_radix`。这 71 处声明已删除，326 处引用直接读取 `orbital_index_radix`。轨道标签计算需要的进制平方继续按公式从中央常量派生。

更新外部调用方时，按当前接口重新编译并链接。CMake 在增量配置时清理已撤下的实现模块和旧 DSPEVX 接口产物，安装时移除以前版本安装的三个旧库模块文件。[公开过程名称对照](../implemented/design-0007-shared-library-procedure-names.md)

第三批验证覆盖 1,058 份 Fortran 源文件。按逻辑语句复扫，旧命名入口、`USE` 重命名和仅用另一个常量名初始化的本地容量同义名均为零。GNU Fortran 15.3 Check（串行、MPI、jj2lsj OpenMP）通过 117/117 项 CTest，纯串行 Release 通过 77/77 项。网格补丁器通过 38 项测试。

安装验证先在临时构建目录和安装前缀放入上一版发布的模块，再重新配置、构建及安装。独立客户端使用正式入口编译链接成功，旧的容量名、过程名、关键字和已撤下的模块在负向编译中均遭拒绝。静态归档及两个 RCI 可执行文件都使用新的 Breit 收集符号。与撤下兼容入口前的 Check 二进制比较，63 次科学输出文件对照逐字节一致。
