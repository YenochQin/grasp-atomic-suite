# 普通程序命名与目录迁移

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

## 应用文件、模块与过程

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
命令行入口沿用下表名称。RCI 的 [Breit 系数收集](../../src/apps/rci/collect_breit_angular_coefficient.f90)
和 [Breit 径向因子](../../src/apps/rci/breit_radial_integral_factors.f90) 使用描述性文件名，
并以链接时选择的 `talk` 和 `cxk` 作为外部符号。
[RMCDHF 来源清单](../reference/rmcdhf_app_symbols_snapshot.json)记录当前路径与历史源码键的对应关系。
应用文件迁移涉及 1,335 个原路径和 1,333 个新路径，其中 650 个是接口文件。
两个未参与 CMake 构建的字节相同的接口副本已合并，另移除了一个使用旧签名的 RCI 接口副本。

命名迁移以 `bb0a2c7` 和当时的工作树为基线。GNU Fortran 15.3.0 的 Release 和 Check 各通过
114 项 CTest。最后校正径向密度输入和单电子占据检查的名称后，两套配置各复验相关 7 项测试。
同一组夹具迁移前后的串行和双进程 MPI 输出各有 16 个文件逐字节一致，RCI 的 Coulomb 和
Breit 混合系数文件也分别逐字节一致。

## 程序目录与可执行文件

日期：2026-11-11。基于 `c47f4d4` 及后续 MPI 描述符输入的工作树修改。
已改造的普通程序采用功能名称，去掉旧 `90` 标记；语言标准由 `f2018-release/check/audit` 构建配置表明。
文件仍使用编译器识别的 `.f90` / `.F90` 扩展名。

| 原 `src/appl` 目录 | 当前 `src/apps` 目录 | 当前可执行文件 |
| --- | --- | --- |
| `gj90` | `gj` | **`gj`、`gj_mpi`**，原名 `gj90`、`gj90_mpi` |
| `rhfs90` | `rhfs` | `rhfs`、`rhfs_mpi` |
| `jj2lsj90` | `jj2lsj` | `jj2lsj_omp` |
| `rangular90` | `rangular` | `rangular` |
| `rangular90_mpi` | `rangular_mpi` | `rangular_mpi` |
| `rwfnestimate90` | `rwfnestimate` | `rwfnestimate` |
| `rci90` | `rci` | `rci` |
| `rci90_mpi` | `rci_mpi` | `rci_mpi` |
| `rbiotransform90` | `rbiotransform` | `rbiotransform` |
| `rbiotransform90_mpi` | `rbiotransform_mpi` | `rbiotransform_mpi` |
| `rtransition90` | `rtransition` | `rtransition` |
| `rtransition90_mpi` | `rtransition_mpi` | `rtransition_mpi` |
| `rtransition90_phase` | `rtransition_phase` | `rtransition_phase` |

重命名阶段曾包含 strict/diag 两个独立求解器，后已按用户确认的否决决定[移除](rejected_rmcdhf_variants.md)。
上表为移除后仍保留的 13 个普通程序目录。

同时将入口文件 `gj90.f90`、`gj90_mpi.f90`、普通 `rci90mpi.f90` 改为
`gj.f90`、`gj_mpi.f90`、`rci_mpi.f90`，并同步程序标识、GJ 帮助信息和 CMake 目标。
其它表中可执行文件原本就没有 `90`，命令不变。

CSFG 继续处于独立可选边界，其目录/程序名称保留。
自旋–角向库已更名为 [`spin_angular`](spin_angular_naming.md)。原 `libmcp90` 已更名为
[`breit_angular_support`](breit_angular_support_migration.md)。
Davidson 求解库已更名为 [`libdavidson_eigensolver`](davidson_eigensolver_migration.md)。
MPI 支持库源码位于 [src/lib/mpi](grasp_mpi_naming.md)。
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
    -DGRASP_BUILD_CSFG=OFF \
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
