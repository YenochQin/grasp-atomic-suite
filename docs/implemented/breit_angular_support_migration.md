# Breit 角系数支持库命名迁移

原 `libmcp90` 的新库名为 `breit_angular_support`，对应目录 `src/lib/breit_angular_support/`、CMake 目标 `breit_angular_support`、
静态库 `libbreit_angular_support.a` 和安装目录 `lib/breit_angular_support/` 下的 Fortran 模块文件。
两个模块分别提供 Breit 角系数缓冲和径向积分系数接口，外部入口的源文件使用 `_external.f90` 后缀。

| 原文件 | 现文件 | Fortran 模块 |
| --- | --- | --- |
| `talk.f90` | `breit_angular_coefficient_buffer.f90` | `talk_i` → `breit_angular_coefficient_buffer` |
| `cxk.f90` | `breit_radial_integral_coefficients.f90` | `cxk_i` → `breit_radial_integral_coefficients` |
| `talk_external.f90` | `breit_angular_coefficient_buffer_external.f90` | 外部 `TALK` 转发入口 |
| `cxk_external.f90` | `breit_radial_integral_coefficients_external.f90` | 外部 `CXK` 转发入口 |

核对日期为 2026-11-20，命名迁移基准为 `2d4c721` 的 `src/lib/libmcp90` 及 RCI 对应实现。基准版本的四份
Fortran 源码共有 33 个声明实体。`TALK` 的 9 个形参和 `CXK` 的 7 个形参分别在过程体和显式接口中
各声明一次，共迁移 32 处声明，另删除 1 个未使用的局部常量。随后合并了库内的默认过程及其接口来源。

## `TALK`

映射适用于 [库模块及接口](../../src/lib/breit_angular_support/breit_angular_coefficient_buffer.f90)、
[外部兼容入口](../../src/lib/breit_angular_support/breit_angular_coefficient_buffer_external.f90)、
[串行 RCI 实现](../../src/apps/rci/talk.f90)和 [MPI RCI 实现](../../src/apps/rci_mpi/talk.f90)。
轨道命名沿用 [SPEAK](../../src/lib/core/mcp_coefficient_output.f90) 的左右态约定。[RCI 矩阵组装](../../src/apps/rci/setham_gg.f90)
按以下顺序读取缓冲区标签，并按类型码选用相应的 Breit 径向积分。

| 原名 | 迁移后 | 含义 |
| --- | --- | --- |
| `ja` | `bra_csf_index` | 左态 CSF 的索引，用于调试行。 |
| `jb` | `ket_csf_index` | 右态 CSF 的索引，用于调试行。 |
| `nu` | `radial_integral_order` | 径向积分的阶数，存入标签第 5 项。 |
| `ia` | `bra_orbital_1` | 左态第一个轨道的索引，存入标签第 1 项。 |
| `ib` | `bra_orbital_2` | 左态第二个轨道的索引，存入标签第 2 项。 |
| `ic` | `ket_orbital_1` | 右态第一个轨道的索引，存入标签第 3 项。 |
| `id` | `ket_orbital_2` | 右态第二个轨道的索引，存入标签第 4 项。 |
| `itype` | `integral_type_code` | 有符号的 Breit 积分类型码，存入标签第 6 项。绝对值选择 `BRINT1` 至 `BRINT6`，负号标记闭芯项。 |
| `coef` | `coefficient` | 与该径向积分相乘的实数角系数。 |

## `CXK`

映射适用于 [库模块及接口](../../src/lib/breit_angular_support/breit_radial_integral_coefficients.f90)、
[外部兼容入口](../../src/lib/breit_angular_support/breit_radial_integral_coefficients_external.f90)、
[串行 RCI 实现](../../src/apps/rci/cxk.f90)和 [MPI RCI 实现](../../src/apps/rci_mpi/cxk.f90)。
形参的含义来自实际计算系数的 RCI 过程，以及构造轨道索引和 `2κ` 数组的 [BREID](../../src/apps/rci/breid.f90)。

| 原名 | 迁移后 | 含义 |
| --- | --- | --- |
| `s` | `radial_integral_coefficients` | Breit 有效相互作用展开中的 12 个径向积分系数槽位。 |
| `is` | `orbital_indices` | 参与相互作用的四个轨道索引。 |
| `kaps` | `twice_orbital_kappa` | 对应四个轨道各自的 `2κ`，其中 `κ` 为 Dirac 角量子数。 |
| `nu` | `radial_integral_order` | 径向积分的阶数。 |
| `k` | `tensor_rank` | 有效相互作用的张量阶数，可与径向积分阶数不同。 |
| `ibr` | `radial_integral_case` | 由轨道索引的相等关系确定的径向积分分类，取值 1 至 4。 |
| `iex` | `direct_exchange_code` | 直接项或交换项的选择码，1 为直接项，2 为交换项。 |

## 局部常量与模块导入

| 原名 | 迁移后 | 用途及依据 |
| --- | --- | --- |
| `talk::key` | 删除 | `orbital_index_radix` 的局部 `PARAMETER` 别名，仅见于声明。`TALK` 分栏存储轨道索引。 |
| `parameter_def::orbital_index_radix` | 移除本地导入 | 此处唯一用途是定义 `key`。 |
| `debug_c::ibug1` | `angular_coefficient_debug_flag` | 本地 `USE` 别名。非零时向 Fortran 单元 99 写一行系数及其积分标签。 |
| `debug_c::cutoff`、`double`、`ibug2` 至 `ibug6`、`ldbpa`、`ldbpg`、`ldbpr` | 移除本地导入 | 库内 `TALK` 没有引用这些共享状态或重导出的 kind。 |
| `buffer_c::double`、`orb_c::double` | 移除本地导入 | 库直接使用标准 kind `real64`。 |
| `buffer_c::angular_buffer_capacity`、`angular_coefficient_count`、`angular_coefficient_labels`、`angular_coefficients` | 沿用原名 | 依次为缓冲区容量、已登记条数、积分标签表和系数表。 |
| `orb_c::orbital_principal_n`、`orbital_symmetry_labels` | 沿用原名 | 调试输出所需的轨道主量子数和对称性标签。 |
| `iso_fortran_env::real64` | 沿用原名 | Fortran 标准 kind 常量，四份源码各导入一次。 |

## 模块接口与链接

| 原接口文件 | 模块内默认实现 | 保留的公开外部入口 |
| --- | --- | --- |
| `talk_I.f90` | `breit_angular_coefficient_buffer::talk_fallback` | `breit_angular_coefficient_buffer::talk`，GNU Fortran 链接符号 `talk_` |
| `cxk_I.f90` | `breit_radial_integral_coefficients::cxk_fallback` | `breit_radial_integral_coefficients::cxk`，GNU Fortran 链接符号 `cxk_` |

默认实现位于各自模块的 `contains` 下。模块用 `procedure(talk_fallback)` 和 `procedure(cxk_fallback)`
从过程定义派生公共外部过程的显式接口。`breit_angular_support` 静态库分别将两个转发入口放在独立归档成员中。
`rci` 和 `rci_mpi` 的同名应用过程提供 `talk_` 和 `cxk_`，链接器因此选用应用实现。通过普通静态库链接而未提供应用实现的调用方，
经外部入口转发至模块中的默认过程。关键字采用上表迁移后的名称，参数顺序、类型、kind、形状和 `intent` 取自同一套接口。

- `spin_angular` 目标内的直接调用方是 `el1`、`el2`、`el31`、`el32`、`el33`、`el41`、`el51`、`el52` 和 `el53`。
  `rci` 和 `rci_mpi` 另各有一份 `breid.f90` 直接调用两者。这 11 份源码共有 177 处调用，均使用位置实参。
  检查的仓库源码中未发现泛型成员、过程实参或额外的 `external` 声明。
- `spin_angular` 依赖 `breit_angular_support`。`rangular`、`rci`、`gj`、`rhfs`、`ris4`、`rtransition` 及其 MPI 版本，
  还有 `rtransition_phase` 直接链接这两个库。`rbiotransform`、`rdensity` 及其 MPI 版本、`rwfnestimate`、
  `jj2lsj` 的内核与转换目标和 `test.f2018_nmtejj` 通过 `spin_angular` 引用这些库。
- 外部 Fortran 工程需要更新 `USE` 模块名、模块搜索路径和链接选项 `-lbreit_angular_support`。公共过程的名称和关键字沿用现有约定，
  使用 `talk_` 或 `cxk_` 外部符号的对象可重新链接新静态库。

库内 `CXK` 仍以原消息停止。RCI 依靠自己的外部过程符号执行 Breit 计算。库内 `TALK` 每次调用都登记系数，
RCI 版 `TALK` 只登记绝对值大于 `eps` 的系数。该差异与调试输出格式均沿用原实现。

## 命名迁移验证

- 使用 GNU Fortran 15.3.0 构建启用 MPI 的 Release 和 Check 配置。Release 全部 110 项 CTest 通过。
  Check 通过 109 项，`rci_mpi_breit` 在 [ZKF](../../src/apps/rci_mpi/zkf.f90#L102) 触发 `SIGFPE`。
  用同一 Check 配置重新链接迁移前的 RCI MPI `TALK` 和 `CXK` 后，复现相同位置的异常。
- 使用现有 7 个 Ni I `J=0+` CSF 的回归夹具，对照迁移前后的 Release 串行和双进程 MPI 程序。
  Breit 启用和关闭的四组对照中，各自的 `case.cm` 逐字节一致。
- 用临时驱动检查库的 Fortran 关键字调用、缓冲区扩容、零系数收集和有符号积分类型码。
  迁移前后的 `TALK` 调试行与占位 `CXK` 的停止消息逐字节一致。

## 模块化验证

- 在 `2d4c721` 的隔离源码树中叠加命名后的 `libmcp90`，用 GNU Fortran 15.3.0 对比模块化前后的 Release 构建。
  启用 MPI 的模块化构建通过全部 110 项 CTest，包含 `rci_mpi_breit`。前后的 `rci`、`rci_mpi` 可执行文件分别逐字节一致。
- 当前工作树的 Release 全量编译通过。`rci_mpi_breit`、`rangular_mpi_descriptors` 和串行、MPI 两项 `f2018_rk_count_rci` 回归全部通过。
- 临时驱动分别以旧版和新版 `.mod` 编译，用相同关键字及过程指针调用。用旧接口编译的对象可重新链接新版归档，
  库外源码也能使用临时安装的模块和静态库。
  两个外部符号可分别由调用方提供，另一过程仍从库中解析。`TALK` 的缓冲区内容和调试输出与旧版一致，`CXK` 仍打印原停止消息。

## 库、文件与模块更名验证

- 使用 GNU Fortran 15.3.0 完成 Release 构建和临时安装。安装生成 `libbreit_angular_support.a` 及两个新模块文件，
  外部驱动使用新模块通过关键字调用。使用旧接口编译的对象可链接新静态库，并产生相同的缓冲区、调试和停止输出。
- 在同一套依赖上对照更名前后的 Release 构建，串行和 MPI RCI 可执行文件分别逐字节一致。
  只叠加本次更名的隔离构建通过全部 110 项 CTest。当前工作树整合其他库名改动后重新构建，
  7 项相关回归全部通过，包含 `rci_mpi_breit` 和 `rangular_mpi_descriptors`。
