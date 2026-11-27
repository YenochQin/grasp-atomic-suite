# 普通程序命名与目录迁移

## 2026-11-21 语义命名审查整改

在 `d9f9787` 的审查基础上，按实际量、编码和计算角色更新了变量、形参和共享状态。
例如 `csf_two_j_plus_one` 表示 RDENSITY 计算使用的 `2J+1`，
`decoded_integral_orbital_indices` 表示 RMCDHF 从积分标签解包的四个轨道索引。
串行和 MPI 的 RMCDHF 轨道势入口现在用同一套形参名称。内存版 MCP 数据区分每项积分的贡献数、
CSF 列、稀疏矩阵位置及角系数。三个跃迁实现使用各积分及贡献项的名称。

| 旧源码相对路径 | 现行源码相对路径 |
| --- | --- |
| `jj2lsj/ls_jj_coefficients_table_2.f90` | [jj2lsj/ls_jj_coefficients_f6.f90](../../src/apps/jj2lsj/ls_jj_coefficients_f6.f90) |
| `jj2lsj/ls_jj_coefficients_table_3.f90` | [jj2lsj/ls_jj_coefficients_f7.f90](../../src/apps/jj2lsj/ls_jj_coefficients_f7.f90) |
| `rangular_mpi/open_temporary_mcp_channel.f90` | [rangular_mpi/open_temporary_mcp_channels.f90](../../src/apps/rangular_mpi/open_temporary_mcp_channels.f90) |
| `ris4/fical.f90` | [ris4/isotope_shift_transition_main.f90](../../src/apps/isotope/isotope_shift_transition_main.f90) |
| `rmcdhf_mem/rmcdhf_main.f90` | [rmcdhf_mem/rmcdhf_mem_main.f90](../../src/apps/rmcdhf_mem/rmcdhf_mem_main.f90) |
| `rmcdhf_mem_mpi/rmcdhf_mpi_main.f90` | [rmcdhf_mem_mpi/rmcdhf_mem_mpi_main.f90](../../src/apps/rmcdhf_mem_mpi/rmcdhf_mem_mpi_main.f90) |
| `rtransition_phase/rtransition_main.f90` | [rtransition_phase/rtransition_phase_main.f90](../../src/apps/rtransition_phase/rtransition_phase_main.f90) |

沿用现行可执行目标和命令，包括 `fical`。`rangular_mpi` 的通道打开过程改为复数
`open_temporary_mcp_channels`，所在模块为 `open_temporary_mcp_channels_interface`。
使用 Fortran 关键字实参的调用方以现行形参声明为准。旧构建清单删除已合并或替换的源文件引用，
CMake 继续作为构建入口。[RMCDHF 来源与改名记录](../reference/rmcdhf_audited_names_snapshot.json)

跃迁工作区的 `jja`、`jjb`、`hb1/2`、`hc1/2`、`hm1/2` 目前只参与分配、扩容和释放。
本轮保留这些指针的旧名，现行计算路径无法区分它们原先的物理角色。

GNU Fortran 15.3.0 和 Open MPI 5.0.11 下，启用串行程序、MPI 和 jj2lsj OpenMP 的
Release 与 Check 配置各通过 117 项 CTest。与 `d9f9787` 的 Release 构建用相同夹具对照，
79 个计算文件或报告输出逐字节一致，覆盖串行和双进程工作流、Coulomb/Breit RCI、RIS4、
跃迁、JJ2LSJ、FICAL 和 RLEVELS。公开 C 头文件的 8 个 ABI 函数名称和类型与基准一致。

## 应用源码目录和变量命名

自 2026-11-20 起，suite 的 20 个应用源码目录位于 `src/apps/`。CMake、测试和网格补丁器的
`atomic-suite` 布局使用该路径。`grasp2018` 布局及历史来源快照中的 `src/appl/` 对应原目录。

应用中的共享变量和形参采用 `src/lib/` 的命名约定。典型用法包括区分 CSF 数量和索引
（`csf_count`、`bra_csf_index`、`ket_csf_index`），区分 MPI rank 和进程数
（`mpi_rank`、`mpi_process_count`），以及区分径向网格和轨道索引
（`radial_index`、`orbital_index`）。实现和显式接口中的同一形参使用同名。使用关键字实参的
Fortran 调用方应以当前接口声明为准。

本次迁移涉及 914 份 Fortran 文件的 3,261 项标识符。在包含当时 I/O 改动的验证副本中，GNU Fortran 15.3.0
的 Release 和 Check 配置均启用串行目标、MPI 和 jj2lsj OpenMP，各通过 114 项 CTest。使用同一组夹具比较迁移前后的
Release 二进制，串行和双进程 MPI 各有 16 个波函数、混合系数和 MCP 文件逐字节一致。
RCI 的 Coulomb 和 Breit 路径各自生成的混合系数文件也逐字节一致。

## 应用文件、模块与过程的命名迁移

应用源码使用小写 snake_case。独立过程文件以过程命名，显式接口文件和模块使用
`<procedure>_interface`。模块文件采用模块名。串行、MPI 和内存版采用同一术语，并用
`_mpi` 区分需要独立符号的 MPI 过程。旧缩写程序入口源码改用 `*_main.f90`。

| 原名称 | 现名称 |
| --- | --- |
| `rscfvu.f90` | [rmcdhf_main.f90](../../src/apps/rmcdhf/rmcdhf_main.f90) |
| `setcof`、`setcof_i` | `assemble_orbital_potential_coefficients`、`assemble_orbital_potential_coefficients_interface` |
| `bndinv`、`bndinv_i` | `invert_banded_matrix`、`invert_banded_matrix_interface` |
| `rmcdhf_mem_c` | `rmcdhf_mcp_memory` |
| `jj2lsj_code` | `jj2lsj_conversion` |
| `mrgcsl` | `merge_transition_csf_bases` |

Fortran 源码调用方需更新 `USE`、过程调用和回调绑定，再重新编译。CMake 的程序目标和
命令行入口沿用下表名称。RCI 在链接时选择应用自己的 [Breit 系数收集过程 `collect_breit_angular_coefficient`](../../src/apps/rci/collect_breit_angular_coefficient.f90)。
串行和 MPI 的 Breit 径向积分系数由共享模块过程
[`calculate_breit_radial_integral_coefficients`](../../src/lib/breit_angular_support/breit_radial_integral_coefficients.f90) 计算。
[RMCDHF 来源清单](../reference/rmcdhf_app_symbols_snapshot.json)记录当前路径与历史源码键的对应关系。
应用文件迁移涉及 1,335 个原路径和 1,333 个新路径，其中 650 个是接口文件。
两个未参与 CMake 构建的字节相同的接口副本已合并，另移除了一个使用旧签名的 RCI 接口副本。

命名迁移以 `bb0a2c7` 和当时的工作树为基线。GNU Fortran 15.3.0 的 Release 和 Check 各通过
114 项 CTest。最后校正径向密度输入和单电子占据检查的名称后，两套配置各复验相关 7 项测试。
同一组夹具迁移前后的串行和双进程 MPI 输出各有 16 个文件逐字节一致，RCI 的 Coulomb 和
Breit 混合系数文件也分别逐字节一致。

## 应用过程模块化

以 `e8a15a9` 为基线，将 647 个应用过程放入原显式接口模块的 `contains` 下。
过程定义提供显式接口，原 `use …, only: …` 调用保留。模块使用实现文件名，
例如 `invert_banded_matrix_interface` 定义在
[invert_banded_matrix.f90](../../src/apps/biorthtransform/invert_banded_matrix.f90)。
RCI 的径向密度比和稀疏矩阵乘积规约分别移到自己的实现文件，MPI 版使用自身的规约实现。
完整的过程、实现文件和原接口对应关系见[迁移清单](../reference/application_procedure_modules.json)。

CMake 目标继续编译实现文件。MCP、CI 和权重测试以同名测试模块提供替身，
并为各测试目标隔离模块输出目录。径向网格容量测试在 790／800／149 的配置下编译
实际的波函数模块。跃迁标签数组的分配过程按用途命名形参，例如 `coefficient_end_offsets`、`right_child_nodes` 和 `sorted_integral_labels`。

最后合并的 13 个过程采用现场调用和实现所需的接口。

| 范围 | 实现提供的接口 |
| --- | --- |
| biorthtransform 和 MPI 版的 `print_matrix_symmetry`、MPI 版的 `write_transformed_radial_wavefunctions` | 用假定大小数组读取打包矩阵，用 128 字符接收两个文件名。 |
| RIS4 的 `fit_electron_density`，rtransition 串行和 phase 版的 `load_initial_lsj_labels` | 输出六项场移电子因子，按假定长度读取 LSJ 文件名。 |
| RCI 串行的 `calculate_vacuum_polarization_integral`、四个 RMCDHF 内存装载过程 | 将会交换的轨道索引定义为 `INTENT(INOUT)`，将用作循环上界的 `ncontr` 定义为 `INTENT(IN)`。 |
| RCI 串行及 MPI 的 `accumulate_closed_shell_breit` | 计算过程保留只读形参。矩阵装配模块通过私有过程适配 `calculate_two_electron_angular_coefficients` 的回调签名。 |

串行和 MPI RCI 直接从 `blas_lapack_interfaces` 导入 `dspevx`。`collect_breit_angular_coefficient` 承担 Breit 系数收集的链接时替换。
已从跃迁目录清理两份无调用者的占据数接口，其实际实现归入 biorthtransform。
suite 的 CMake 流程使用新模块。工作区其他源码仓库未发现对已迁移外部符号的调用。
如有构建树之外的调用方直接链接这些旧外部符号，应改用相应模块并重新编译和链接。

GNU Fortran 15.3.0 的干净 Release 和 Check 构建均启用串行、MPI 和 jj2lsj OpenMP。
最后一批合并后，两套配置各通过 116 项 CTest。首批合并逐项比较了 742 个按生产目标展开的编译接口，
并单独编译和核对了 3 个构建清单外的过程。串行和双进程 MPI 夹具各有 16 个波函数、
混合系数及 MCP 文件逐字节一致。RCI 的 Coulomb 和 Breit 混合系数文件也分别逐字节一致。
专项对照还确认了串行 RCI 的真空极化混合系数和 RIS4 的六项场移因子逐字节一致。

## 程序目录与可执行文件

日期：2026-11-11。基于 `c47f4d4` 及后续 MPI 描述符输入的工作树修改。
已改造的普通程序采用功能名称，去掉旧 `90` 标记；语言标准由 `f2018-release/check/audit` 构建配置表明。
文件仍使用编译器识别的 `.f90` / `.F90` 扩展名。

| 原 `src/appl` 目录 | 当前 `src/apps` 目录 | 当前可执行文件 |
| --- | --- | --- |
| `gj90` | `gj` | **`gj`、`gj_mpi`**，原名 `gj90`、`gj90_mpi` |
| `rhfs90` | `hyperfine` | `hyperfine`、`hyperfine_mpi`，此前为 `rhfs`、`rhfs_mpi` |
| `ris4` | `isotope` | `isotope`、`isotope_mpi`，此前为 `ris4`、`ris4_mpi` |
| `jj2lsj90` | `jj2lsj` | `jj2lsj_omp` |
| `rangular90` | `rangular` | `rangular` |
| `rangular90_mpi` | `rangular_mpi` | `rangular_mpi` |
| `rwfnestimate90` | `rwfnestimate` | `rwfnestimate` |
| `rci90` | `rci` | `rci` |
| `rci90_mpi` | `rci_mpi` | `rci_mpi` |
| `rbiotransform90` | `biorthtransform` | `biorthtransform` |
| `rbiotransform90_mpi` | `biorthtransform_mpi` | `biorthtransform_mpi` |
| `rtransition90` | `rtransition` | `rtransition` |
| `rtransition90_mpi` | `rtransition_mpi` | `rtransition_mpi` |
| `rtransition90_phase` | `rtransition_phase` | `rtransition_phase` |

重命名阶段曾包含 strict/diag 两个独立求解器，后已按用户确认的否决决定[移除](adr-0002-retire-strict-and-diagnostic-rmcdhf.md)。
双正交变换入口经 `rbiotransform`、`rbiotransform_mpi` 迁移至表中的现名。

同时将入口文件 `gj90.f90`、`gj90_mpi.f90`、普通 `rci90mpi.f90` 改为
`gj.f90`、`gj_mpi.f90`、`rci_mpi.f90`，并同步程序标识、GJ 帮助信息和 CMake 目标。
其它表中可执行文件原本就没有 `90`，命令不变。

自旋–角向库已更名为 [`spin_angular`](design-0012-spin-angular-library-naming.md)。原 `libmcp90` 已更名为
[`breit_angular_support`](design-0015-breit-angular-support-migration.md)。
Davidson 求解库已更名为 [`libdavidson_eigensolver`](design-0014-davidson-eigensolver-migration.md)。
MPI 支持库源码位于 [src/lib/mpi](design-0010-mpi-library-naming.md)。
网格补丁器的 `atomic-suite` 布局使用新路径，`grasp2018` 布局保留上游原路径。
来源快照 JSON、历史报告中的原始命令和已记录身份属于当时版本；指向现行源码的文档链接已更新。

F2018 验收清单和警告基线同步迁移文件路径、测试/目标名称及诊断身份散列；
各配置的编译数量、已有警告类别和数量预算不因命名迁移而改变。
这次命名调整不改变数值算法或宣称完成尚未验收的 F2018/物理验证。

## 更新后的构建

建议使用新构建目录，避免原构建目录或安装目录中残留旧 `gj90` 文件。
下面的普通版构建包含串行、MPI 及已有 OpenMP 目标：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
cmake -S . -B build-f2018-programs \
    -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
    -DCMAKE_BUILD_TYPE=Release -DGRASP_FORTRAN_PROFILE=f2018-release \
    -DGRASP_BUILD_SERIAL_APPS=ON -DGRASP_ENABLE_MPI=ON \
    -DGRASP_JJ2LSJ_OPENMP=ON
cmake --build build-f2018-programs -j4
```

调用 Landé 因子程序的用户脚本应使用 `gj` / `gj_mpi`；本次只更新 suite 仓库，
工作区其它工具仓库的外部配置由各自维护。本次未创建旧名称别名或删除既有安装文件。

## 验证记录

重命名提交 `212bcc6` 时，GNU Fortran 15.3.0、Open MPI 5.0.11 下的普通版配置在新目录中完整构建通过，
当时的 79 项 CTest 全部通过，F2018 目标/测试清单检查通过。后续退役后的验证见移除说明。
网格补丁回归同时覆盖上游和 suite 两种布局；现行源码文档链接已核对。
命名调整中的 Fortran 内容变化仅涉及三个入口文件的程序标识/GJ 名称与提示，
其余移动的 Fortran 源文件与移动前内容一致。
本轮未重跑全部构建配置矩阵，也未重复性能测量或进行新的物理收敛实验。
本地日志保存在忽略目录 `build-program-rename-20261111/`。
