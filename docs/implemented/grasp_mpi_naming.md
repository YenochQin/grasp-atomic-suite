# GRASP MPI 支持库命名迁移

核对日期为 2026-11-20，源码基准为 `2d4c721` 上已经完成变量命名迁移的工作树。
源码目录从 `src/lib/mpi90` 改为 `src/lib/mpi`，CMake 目标从 `mpiu90` 改为 `grasp_mpi`，
静态库安装为 `libgrasp_mpi.a`。文件名与模块名逐一对应，以 `mpi_` 开头，后缀表示各模块职责。

## 文件与模块对照

| 原文件 | 原模块 | 现模块与源码 | 职责 |
| --- | --- | --- | --- |
| `mpi_C.f90` | `mpi_c` | [`mpi_state`](../../src/lib/mpi/mpi_state.f90) | 共享 MPI 状态和 MPI 接口重导出。 |
| `mpi_runtime.f90` | `mpi_runtime_i` | [`mpi_runtime`](../../src/lib/mpi/mpi_runtime.f90) | MPI 启动、结束与连续字节通信类型。 |
| `mpiu.f90` | `mpiu_i` | [`mpi_utilities`](../../src/lib/mpi/mpi_utilities.f90) | 程序启停、计时、消息汇总与全局归约。 |
| `cpath.f90` | `cpath_i` | [`mpi_work_paths`](../../src/lib/mpi/mpi_work_paths.f90) | 设置各 rank 的临时目录与 rank 0 的串行 I/O 路径。 |
| `sys_chdir.f90` | `sys_chdir_i` | [`mpi_change_directory`](../../src/lib/mpi/mpi_change_directory.f90) | 切换进程工作目录。 |
| `sys_getwd.f90` | `sys_getwd_i` | [`mpi_current_directory`](../../src/lib/mpi/mpi_current_directory.f90) | 获取进程工作目录。 |
| `sys_mkdir.f90` | `sys_mkdir_i` | [`mpi_make_directory`](../../src/lib/mpi/mpi_make_directory.f90) | 创建目录。 |
| `input_io_mpi.f90` | `input_io_mpi` | [`mpi_io_errors`](../../src/lib/mpi/mpi_io_errors.f90) | 报告输入错误并中止 MPI 通信域。 |
| `csf_basis_mpi.f90` | `csf_basis_mpi` | [`mpi_csf_basis_input`](../../src/lib/mpi/mpi_csf_basis_input.f90) | 加载并广播完整 CSF 基组。 |
| `csf_block_mpi_input.f90` | `csf_block_mpi_input` | [`mpi_csf_block_input`](../../src/lib/mpi/mpi_csf_block_input.f90) | 为 RCI 按对称性块加载 CSF。 |
| `csf_descriptor_collectives.f90` | `csf_descriptor_collectives` | [`mpi_csf_collectives`](../../src/lib/mpi/mpi_csf_collectives.f90) | 同步描述符状态并按批广播 CSF 数组。 |
| `rmcdhf_descriptor_mpi.f90` | `rmcdhf_descriptor_mpi` | [`mpi_rmcdhf_csf_input`](../../src/lib/mpi/mpi_rmcdhf_csf_input.f90) | 加载 RMCDHF 的描述符并核对对称性块。 |
| `setisompi.f90` | `setisompi_i` | [`mpi_isotope_input`](../../src/lib/mpi/mpi_isotope_input.f90) | 加载并广播核质量、核模型与核矩参数。 |
| `setrwfmpi.f90` | `setrwfmpi_i` | [`mpi_radial_wavefunction_file`](../../src/lib/mpi/mpi_radial_wavefunction_file.f90) | 打开径向波函数文件并校验文件头。 |
| `lodrwfmpi.f90` | `lodrwfmpi_i` | [`mpi_radial_wavefunction_reader`](../../src/lib/mpi/mpi_radial_wavefunction_reader.f90) | 读入、广播和插值径向轨道。 |
| `npol_init_mpi.f90` | `npol_mpi` | [`mpi_nuclear_polarization`](../../src/lib/mpi/mpi_nuclear_polarization.f90) | 同步核极化配置、势和初始化状态。 |
| `iniestmpi.f90` | `iniestmpi_i` | [`mpi_initial_eigenvectors`](../../src/lib/mpi/mpi_initial_eigenvectors.f90) | 用稀疏子矩阵构造 Davidson 初始向量。 |
| `spicmvmpi.f90` | `spicmvmpi_i` | [`mpi_sparse_matrix_product`](../../src/lib/mpi/mpi_sparse_matrix_product.f90) | 计算各 rank 的稀疏哈密顿量乘积并归约。 |

## 调用与构建

- Fortran 调用方按对照表更新 `use` 语句，并用新模块重新编译。过程名称和形参沿用现有接口。
- CMake 调用方链接 `grasp_mpi`。模块文件安装到 CMake 配置前缀下的 `lib/grasp_mpi/`，静态库位于 `lib/`。
- 稀疏矩阵回归测试改名为 `grasp_mpi_local_sparse`。变量、常量及关键字参数见 [变量迁移表](mpi90_variable_migration.md)。

## 验证

- GNU Fortran 15.3.0 与 Open MPI 的 Release 构建通过全部 110 项 CTest。关闭 MPI 时，CMake 目标清单符合串行配置。
- 临时安装生成了 `libgrasp_mpi.a` 和 18 个 `mpi_*.mod` 文件。用安装产物编译的独立程序通过了双进程归约测试。
- RMCDHF 来源检查通过，覆盖 480 个文件。
