# libgrasp_mpi 变量与常量命名迁移

核对日期为 2026-11-20。声明盘点基准为 `2d4c721` 的 `src/lib/mpi90`。实施时纳入工作区已进行的
`double` 导入清理。
每行对应一个作用域中的声明，包含模块变量、局部变量、过程形参、函数结果和命名常量。
Fortran 标识符不区分大小写，表中采用小写。现行文件与模块见 [MPI 库命名迁移](grasp_mpi_naming.md)。

## 命名约定

- `rank` 是通信域内的进程编号，`process_count` 是进程数。`mpi_return_code` 用于接收 MPI 例程返回码，
  I/O、描述符读入、目录操作及 LAPACK 状态按来源命名。`root_` 表示 rank 0 提供的数据。
- `*_index` 是数组索引，`*_count` 是元素数，`*_length` 是有效字符数，`*_capacity` 是缓冲区容量。
  `subshell_occupations`、`subshell_quantum_codes`、`coupling_quantum_codes` 对应三类 CSF 数据。
  `csf_batch_limit` 是每次广播的 CSF 列数上限。不同入口的 `basis_counts` 槽位顺序见对应表。
- 径向波函数和核极化变量沿用 [libcore 命名](core_variable_migration.md)，Davidson 的初始向量和
  工作区沿用 [Davidson 求解库命名](davidson_eigensolver_migration.md)。`—` 表示移除无引用的局部声明。
- Fortran 与 MPI 标准符号沿用公开名称，例如 `real64`、`int8`、`int64`、`error_unit`、
  `MPI_COMM_WORLD`、`MPI_INTEGER1`、`MPI_STATUS_SIZE`。`mpi_state` 继续重导出 MPI 接口。
  历史公开变量 `ierrtotal` 已随旧命名接口移除。
- [CMakeLists.txt](../../src/lib/mpi/CMakeLists.txt) 无本地命名变量或常量声明。`MPI_Fortran_INCLUDE_PATH`、
  `MPI_Fortran_LIBRARIES`、`MPI_Fortran_COMPILE_FLAGS`、`MPI_Fortran_LINK_FLAGS` 由 FindMPI 提供，
  `CMAKE_Fortran_COMPILER_ID` 由 CMake 提供。
- 过程名、位置参数顺序、声明类型、数组布局、MPI 消息和文件格式沿用原约定。若以 Fortran 关键字调用，
  请按表使用新形参名。仓库内使用旧关键字 `column_limit=` 的
  [描述符集合通信测试](../../test/descriptor_collectives.f90)已迁移到 `csf_batch_limit=`。

18 份 Fortran 源码共登记 **294** 个声明实体，重命名 **251** 个，沿用 **25** 个，移除 **18** 个未使用的局部变量。

## [mpi_work_paths.f90](../../src/lib/mpi/mpi_work_paths.f90)

### `cpath`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `startdir` | `startup_directory` | 各进程进入程序时的工作目录。 |
| `permdir` | `serial_io_directory` | rank 0 用于串行 I/O 的目录。 |
| `tmpdir` | `rank_scratch_root` | 返回当前 rank 的临时根目录，其下的编号子目录是实际工作目录。 |
| `ierr_all_nodes` | `directory_status_sum` | 各 rank 目录操作状态的求和结果。 |
| `ierr2m` | `mpi_return_code` | 通信调用返回码。 |
| `lendisk0` | `scratch_path_capacity` | 路径通信缓冲区的固定字符数，128。 |
| `lenidstring` | `rank_directory_name_width` | rank 子目录名宽度，3，对应 I3.3 格式。 |
| `disk` | `scratch_path` | 读入或收发临时根目录路径的固定长度缓冲区。 |
| `idstring` | `rank_directory_name` | 以零补齐的 rank 子目录名。 |
| `lendisk` | `scratch_path_length` | 临时目录路径的有效字符数。 |
| `i` | `destination_rank` | 接收工作路径的 rank。 |
| `len_cwd` | `startup_directory_length` | 原工作目录的有效字符数。 |
| `lstring` | — | 仅声明，无引用。 |
| `asleep` | `retry_delay_value` | 重试间隔计算的暂存值，未传给 sleep。 |
| `itry` | `retry_attempt` | 创建目录或切换目录时的重试序号。 |
| `iteery` | `delay_iteration` | 随 rank 变化的重试间隔循环计数器。 |
| `ydisks` | `disks_file_exists` | 启动目录是否有路径配置文件 disks。 |
| `mpi_tmp` | `scratch_environment_path` | MPI_TMP 的值或回落的 /scratch/USER 路径。 |
| `env_user` | `user_name` | USER 环境变量读出的用户名。 |

### `getenv_allocating`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `getenv_allocating` | `is_defined` | 显式逻辑函数结果，表示是否取得环境变量。 |
| `variable_name` | `variable_name` | 环境变量名。 |
| `value` | `value` | 按需分配的环境变量值。 |
| `length` | `value_length` | 查询得到的字符数。 |
| `status` | `environment_status` | GET_ENVIRONMENT_VARIABLE 返回状态。 |
| `test` | `probe_value` | 用于第一次查询长度的单字符缓冲。 |

## [mpi_csf_basis_input.f90](../../src/lib/mpi/mpi_csf_basis_input.f90)

### `load_csf_basis_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `name` | `descriptor_path` | CSF 描述符路径。 |
| `ncore` | `core_subshell_count` | 闭核子壳层数。 |
| `status` | `descriptor_status` | 描述符读入或工作进程分配状态。 |
| `ierr` | `mpi_return_code` | MPI 广播返回码。 |
| `info` | `basis_counts` | 按 CSF 数、轨道数、电子数、块数、闭核子壳层数排列的广播载荷。 |

## [mpi_csf_block_input.f90](../../src/lib/mpi/mpi_csf_block_input.f90)

### `open_csf_blocks_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `name` | `descriptor_path` | RCI 的 CSF 描述符路径。 |
| `ncore` | `core_subshell_count` | 闭核子壳层数。 |
| `labels` | `block_labels` | 对称性块标签数组。 |
| `status` | `descriptor_status` | 读入描述符块表或分配计数数组的状态。 |
| `ierr` | `mpi_return_code` | MPI 广播返回码。 |
| `info` | `basis_counts` | 按 CSF 数、轨道数、闭核子壳层数、电子数、块数排列的广播载荷。 |

### `load_csf_block_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `block_index` | `block_index` | 本次需要加载的对称性块序号。 |
| `status` | `descriptor_status` | rank 0 读入指定块的状态。 |

## [mpi_csf_collectives.f90](../../src/lib/mpi/mpi_csf_collectives.f90)

### `check_descriptor_status`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `status` | `status` | 调用者提供的局部返回码，非零表示失败。 |
| `failed` | `rank_failed` | 当前 rank 的 0 或 1 失败标志。 |
| `combined` | `any_rank_failed` | 跨 rank 取最大值得到的失败标志。 |
| `ierr` | `mpi_return_code` | MPI_Allreduce 返回码。 |

### `allocate_descriptor_arrays`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `rows` | `csf_count` | 所分配数组的 CSF 数，来自描述符记录数。 |
| `occupations` | `subshell_occupations` | 每个 CSF 的子壳层占据数。 |
| `quantum` | `subshell_quantum_codes` | 每个 CSF 的子壳层量子数编码，每个轨道三项。 |
| `couplings` | `coupling_quantum_codes` | 每个 CSF 的耦合量子数编码。 |
| `status` | `allocation_status` | 当前 rank 的数组分配状态。 |

### `broadcast_descriptor_arrays`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `occupations` | `subshell_occupations` | 每个 CSF 的子壳层占据数。 |
| `quantum` | `subshell_quantum_codes` | 可选的子壳层量子数编码。 |
| `couplings` | `coupling_quantum_codes` | 可选的耦合量子数编码。 |
| `column_limit` | `csf_batch_limit` | 可选的每次广播列数上限，一列对应一个 CSF。 |
| `dimensions` | `local_array_layout` | 本地轨道数、CSF 数及两个可选数组的存在标志。 |
| `root_dimensions` | `root_array_layout` | rank 0 的布局。校验后第一个元素复用为分批列数。 |
| `first` | `first_csf_index` | 本次广播的首列索引。 |
| `last` | `last_csf_index` | 本次广播的末列索引。 |
| `columns` | `csfs_per_batch` | 每次 MPI 广播最多传输的 CSF 数。 |
| `ierr` | `mpi_return_code` | MPI 广播返回码。 |
| `status` | `layout_status` | 本 rank 的维度与可选数组形状检查状态。 |
| `width` | `max_elements_per_csf` | 单个广播数组中每个 CSF 最多贡献的元素数，用于限制 MPI count。 |

## [mpi_initial_eigenvectors.f90](../../src/lib/mpi/mpi_initial_eigenvectors.f90)

### `iniestmpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `nmax` | `submatrix_order_limit` | 用于生成初始向量的子矩阵阶数上限。 |
| `ncf` | `csf_count` | 全哈密顿量阶数及每个全空间向量的长度。 |
| `niv` | `initial_vector_count` | 请求求出的初始近似本征向量数。 |
| `basis` | `basis_vectors` | 全空间初始向量及尾部的对应本征值。 |
| `hmx` | `hamiltonian_values` | 当前 rank 存储的稀疏矩阵元。 |
| `jcol` | `stored_segment_end_offsets` | 按列索引的稀疏段结束位置。 |
| `irow` | `matrix_partner_indices` | 稀疏矩阵元的配对行索引。 |
| `ns` | `submatrix_order` | 实际选取的子矩阵阶数。 |
| `joffspar` | `sparse_segment_offset` | 已处理的本地稀疏存储末尾偏移。 |
| `j` | `column_index` | 遍历原矩阵列，回写时索引本征向量列。 |
| `joffnorm` | `packed_column_offset` | 当前列在打包上三角矩阵中的起始偏移。 |
| `ir` | `stored_element_index` | 本地稀疏存储中的元素索引。 |
| `nfound` | `computed_pair_count` | LAPACK 报告实际算出的本征对数。 |
| `info` | `lapack_status` | DSPEVX 的 INFO 返回码。 |
| `iwork` | `lapack_integer_workspace` | DSPEVX 整数工作区。 |
| `ifail` | `failed_vector_indices` | DSPEVX 报告未收敛的向量索引。 |
| `ap` | `packed_submatrix` | 汇总各 rank 贡献后供 DSPEVX 求解的打包实对称子矩阵。 |
| `eigval` | `initial_eigenvalues` | 子矩阵的低端本征值。 |
| `vec` | `submatrix_eigenvectors` | 子矩阵的本征向量。 |
| `work` | `lapack_real_workspace` | DSPEVX 实数工作区。 |
| `abstol` | `lapack_abs_tolerance` | DSPEVX 的绝对容差。 |

## [mpi_io_errors.f90](../../src/lib/mpi/mpi_io_errors.f90)

### `fail_input_io_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `status` | `io_status` | 调用者提供的 I/O 失败状态。 |
| `unit` | `io_unit` | 发生失败的 Fortran 逻辑单元。 |
| `context` | `context` | 读入错误发生的上下文。 |

### `abort_input_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `status` | `mpi_return_code` | MPI 状态查询或中止调用的返回码。 |
| `initialized` | `mpi_is_initialized` | 是否已调用 MPI 初始化。 |
| `finalized` | `mpi_is_finalized` | 是否已结束 MPI。 |

## [mpi_radial_wavefunction_reader.f90](../../src/lib/mpi/mpi_radial_wavefunction_reader.f90)

### `lodrwfmpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `ierror` | `load_status` | 返回给调用者的轨道加载状态。 |
| `j` | `orbital_index` | 套件中目标轨道的索引。 |
| `i` | `radial_index` | 径向采样点索引。 |
| `k` | `abs_orbital_kappa` | 轨道相对论量子数 κ 的绝对值。 |
| `nwin` | `loaded_orbital_count` | 已匹配并插值的目标轨道数。 |
| `ios` | `read_status` | rank 0 读入轨道记录头的 IOSTAT。 |
| `npy` | `file_principal_n` | 轨道文件记录中的主量子数。 |
| `naky` | `file_kappa` | 轨道文件记录中的相对论量子数 κ。 |
| `my` | `file_grid_size` | 记录中的径向网格点数。 |
| `con` | `squared_charge_to_light_speed_ratio` | (Z/c)²，用于计算点核的轨道近原点指数。 |
| `fkk` | `kappa_squared` | κ² 的实数值。 |
| `ey` | `file_energy_parameter` | 轨道记录中的能量参数。 |
| `pzy` | `file_origin_coefficient` | 轨道记录中的近原点系数。 |
| `dnorm` | `original_norm_squared` | 插值后、归一化前的径向范数平方，由 INTRPQ 返回。 |
| `pa` | `file_large_component` | 读入的大分量数组。 |
| `qa` | `file_small_component` | 读入的小分量数组。 |
| `ra` | `file_radial_grid` | 读入的径向网格。 |

## [mpi_state.f90](../../src/lib/mpi/mpi_state.f90)

### 模块变量

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `processor_name` | `processor_name` | 本 rank 的处理器名称。 |
| `mpi_rank` | `mpi_rank` | 当前 rank 编号。 |
| `mpi_process_count` | `mpi_process_count` | MPI_COMM_WORLD 中的进程数。 |
| `processor_name_length` | `processor_name_length` | 处理器名称的有效字符数。 |
| `mpi_error_code` | `mpi_error_code` | 套件程序共用的整型状态槽。 |
| `ierrtotal` | — | 已移除，仅用于旧接口的公开变量。 |
| `mpi_status_buffer` | `mpi_status_buffer` | MPI_Recv 的整数状态缓冲区。 |

## [mpi_runtime.f90](../../src/lib/mpi/mpi_runtime.f90)

### `startmpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `myid` | `rank` | 当前进程的 rank。 |
| `nprocs` | `process_count` | MPI_COMM_WORLD 的进程数。 |
| `host` | `processor_name` | MPI 返回的处理器名称。 |
| `lenhost` | `processor_name_length` | 处理器名称的有效字符数。 |
| `ierr` | `mpi_return_code` | MPI 初始化和查询调用的返回码。 |

### `mpix_bytes`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `byte_count` | MPI 派生类型中的连续字节数。 |
| `newtype` | `byte_datatype` | 已提交的 MPI 类型句柄。 |
| `ierr` | `mpi_return_code` | 返回给调用方的 MPI_Type_contiguous 和 MPI_Type_commit 状态之和。 |
| `ier0` | `create_status` | MPI_Type_contiguous 的返回码。 |

### `stopmpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `what` | `message` | rank 0 最终打印的消息。 |
| `myid` | `rank` | 当前进程的 rank。 |
| `ierr` | `mpi_return_code` | MPI 屏障和结束调用的返回码。 |

## [mpi_utilities.f90](../../src/lib/mpi/mpi_utilities.f90)

### `mpix_startup`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `myid` | `rank` | 当前进程的 rank。 |
| `nprocs` | `process_count` | MPI_COMM_WORLD 的进程数。 |
| `host` | `processor_name` | MPI 返回的处理器名称。 |
| `lenhost` | `processor_name_length` | 处理器名称的有效字符数。 |
| `ncount1` | `start_ticks` | 开始时间的 SYSTEM_CLOCK 计数。 |
| `progname` | `program_name` | 启动信息中的程序名。 |
| `description` | `description` | 启动信息中的程序说明。 |
| `infiledesc` | `input_file_description` | 控制台显示的输入文件说明。 |
| `outfiledesc` | `output_file_description` | 控制台显示的输出文件说明。 |
| `quiet_workers` | `quiet_workers` | 是否将工作 rank 的标准输出重定向到 /dev/null。 |
| `ncount_rate` | `ticks_per_second` | SYSTEM_CLOCK 返回的每秒计数。 |
| `ncount_max` | `clock_count_limit` | SYSTEM_CLOCK 返回的最大计数。 |
| `iouterr` | `output_open_status` | 重定向标准输出时的 IOSTAT。 |
| `output_message` | `output_message` | 重定向失败时的 IOMSG 缓冲。 |

### `startmpi2`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `myid` | `rank` | 当前进程的 rank。 |
| `nprocs` | `process_count` | MPI_COMM_WORLD 的进程数。 |
| `host` | `processor_name` | MPI 返回的处理器名称。 |
| `lenhost` | `processor_name_length` | 处理器名称的有效字符数。 |
| `ncount1` | `start_ticks` | 工作目录就绪后的 SYSTEM_CLOCK 起始计数。 |
| `startdir` | `startup_directory` | 当前 rank 启动时的工作目录。 |
| `permdir` | `serial_io_directory` | rank 0 用于串行 I/O 的目录。 |
| `tmpdir` | `rank_scratch_root` | 各 rank 的临时根目录。 |
| `progname` | `program_name` | 启动信息中的程序名。 |
| `idstring` | `rank_label` | 以零补齐的三字符 rank 标签。 |
| `chdate` | `date_text` | DATE_AND_TIME 的日期字符串。 |
| `chtime` | `time_text` | DATE_AND_TIME 的时间字符串。 |
| `chzone` | `zone_text` | DATE_AND_TIME 的 UTC 偏移字符串。 |
| `nymduhmsm` | `date_time_values` | DATE_AND_TIME 的整数分量数组。 |
| `lenstart` | `startup_directory_length` | 原工作目录的有效字符数。 |
| `lenperm` | `serial_io_directory_length` | 串行 I/O 目录的有效字符数。 |
| `lentmp` | `rank_scratch_root_length` | rank 临时根目录的有效字符数。 |
| `ncount_rate` | `ticks_per_second` | SYSTEM_CLOCK 返回的每秒计数。 |
| `ncount_max` | `clock_count_limit` | SYSTEM_CLOCK 返回的最大计数。 |
| `msg` | `message` | 发往 rank 0 的启动信息缓冲。 |

### `mpix_shutdown`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `myid` | `rank` | 当前进程的 rank。 |
| `ncount1` | `start_ticks` | 启动时的 SYSTEM_CLOCK 计数。 |
| `progname` | `program_name` | 完成信息中的程序名。 |
| `ncount2` | `elapsed_ticks` | 先接收结束计数，再原位减去起始计数。 |
| `ncount_rate` | `ticks_per_second` | SYSTEM_CLOCK 返回的每秒计数。 |
| `ncount_max` | `clock_count_limit` | SYSTEM_CLOCK 返回的最大计数。 |
| `nseconds` | `elapsed_seconds` | 以整秒表示的运行时间。 |
| `ierr` | `mpi_return_code` | MPI_Finalize 返回码。 |
| `chdate` | `date_text` | 结束日期字符串。 |
| `chtime` | `time_text` | 结束时间字符串。 |
| `chzone` | `zone_text` | 结束时的 UTC 偏移字符串。 |
| `nymduhmsm` | `date_time_values` | DATE_AND_TIME 的整数分量数组。 |

### `stopmpi2`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `myid` | `rank` | 当前进程的 rank。 |
| `nprocs` | `process_count` | MPI_COMM_WORLD 的进程数。 |
| `host` | `processor_name` | MPI 返回的处理器名称。 |
| `lenhost` | `processor_name_length` | 处理器名称的有效字符数。 |
| `ncount1` | `start_ticks` | 启动时的 SYSTEM_CLOCK 计数。 |
| `progname` | `program_name` | 完成信息中的程序名。 |
| `ncount2` | `elapsed_ticks` | 先接收结束计数，再原位减去起始计数。 |
| `ncount_rate` | `ticks_per_second` | SYSTEM_CLOCK 返回的每秒计数。 |
| `ncount_max` | `clock_count_limit` | SYSTEM_CLOCK 返回的最大计数。 |
| `nseconds` | `elapsed_seconds` | 以整秒表示的运行时间。 |
| `chdate` | `date_text` | 结束日期字符串。 |
| `chtime` | `time_text` | 结束时间字符串。 |
| `chzone` | `zone_text` | 结束时的 UTC 偏移字符串。 |
| `nymduhmsm` | `date_time_values` | DATE_AND_TIME 的整数分量数组。 |
| `str2nds` | `elapsed_seconds_text` | 格式化后的整秒数字。 |
| `msg` | `message` | 发往 rank 0 的结束信息缓冲。 |
| `ierr` | `mpi_return_code` | 结束前 MPI 屏障的返回码。 |

### `mpix_printmsg`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `msg` | `message` | 本 rank 的待打印信息，rank 0 复用此缓冲接收工作 rank 的消息。 |
| `myid` | `rank` | 当前进程的 rank。 |
| `nprocs` | `process_count` | MPI_COMM_WORLD 的进程数，也用于构造消息标签。 |
| `inid` | `source_rank` | rank 0 接收消息时遍历的发送方 rank。 |
| `istat` | `receive_status` | MPI_Recv 写入的整数状态数组。 |
| `ierr` | `mpi_return_code` | 消息发送或接收调用的返回码。 |
| `msglength` | `message_length` | 消息的有效字符数。 |

### `mpix_chkpt`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `myid` | `rank` | 当前进程的 rank。 |
| `what` | `checkpoint_message` | 显示给 rank 0 用户的检查点提示。 |
| `ierr` | `mpi_return_code` | 检查点 MPI 屏障的返回码。 |
| `chtmp` | `confirmation_character` | rank 0 从标准输入读取的单字符。 |

### `gdsummpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `x` | `values` | 局部实数输入，返回时各 rank 都得到逐元素求和结果。 |
| `n` | `value_count` | 待归约的元素数。 |
| `ierr` | `mpi_return_code` | MPI_Allreduce 返回码。 |
| `y` | `summed_values` | 接收全局实数求和结果的缓冲。 |

### `gdrsummpi_root`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `x` | `values` | 局部实数输入，rank 0 原位接收逐元素求和结果。 |
| `n` | `value_count` | 待归约的元素数。 |
| `ierr` | `mpi_return_code` | MPI_Reduce 返回码。 |
| `y` | `summed_values` | rank 0 接收实数求和结果的缓冲。 |

### `gdmaxmpi_root`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `x` | `values` | 局部实数输入，rank 0 原位接收逐元素最大值。 |
| `n` | `value_count` | 待归约的元素数。 |
| `ierr` | `mpi_return_code` | MPI_Reduce 返回码。 |
| `y` | `maximum_values` | rank 0 接收实数最大值的缓冲。 |

### `gisummpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `ix` | `values` | 局部整数输入，返回时各 rank 都得到逐元素求和结果。 |
| `n` | `value_count` | 待归约的元素数。 |
| `ierr` | `mpi_return_code` | MPI_Allreduce 返回码。 |
| `i` | — | 仅声明，无引用。 |
| `iy` | `summed_values` | 接收全局整数求和结果的缓冲。 |

## [mpi_nuclear_polarization.f90](../../src/lib/mpi/mpi_nuclear_polarization.f90)

### `npol_init_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `path` | `config_path` | 核极化配置文件路径。 |
| `z` | `nuclear_charge` | 用于核对配置的核电荷数。 |
| `n` | `grid_point_count` | 当前启用的径向网格点数。 |
| `r` | `radial_grid` | 径向网格坐标。 |
| `zz` | `negative_radius_times_potential` | -rV 数组，初始化成功后包含核极化修正。 |
| `comm` | `communicator` | 参加核极化初始化的 MPI 通信域。 |
| `status` | `status` | 在各 rank 同步后的初始化状态。 |
| `message` | `message` | 返回给调用者的错误信息。 |
| `cfg` | `config` | 从 rank 0 读入并分发的核极化配置。 |
| `a` | `dipole_potential` | E1 电偶极极化势。 |
| `b` | `quadrupole_potential` | E2 电四极极化势。 |
| `v` | `polarization_potential` | 总核极化势。 |
| `candidate` | `composed_potential` | 核极化修正后的 -rV，检查通过后写回。 |
| `rank` | `rank` | 当前进程在指定通信域中的编号。 |
| `ierr` | `mpi_return_code` | MPI 查询及广播返回码。 |
| `root_n` | `root_grid_point_count` | rank 0 提供的网格点数，用于检查各 rank 是否一致。 |
| `detail` | `status_message` | 固定长度的通信错误信息缓冲。 |

### `sync_status`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `comm` | `communicator` | 参加错误同步的 MPI 通信域。 |
| `rank` | `rank` | 当前进程在指定通信域中的编号。 |
| `status` | `status` | 进入时是局部状态，返回时是共同状态。 |
| `message` | `message` | 出错时接收首个失败 rank 的错误信息。 |
| `failed_rank` | `failed_rank` | 当前 rank 出错时取其编号，否则取哨兵值。 |
| `first_rank` | `first_failed_rank` | 取最小值后得到的第一个失败 rank，或表示无错误的哨兵。 |
| `ierr` | `mpi_return_code` | MPI 归约或广播返回码。 |

### `broadcast_config`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `cfg` | `config` | 待广播的核极化配置。 |
| `comm` | `communicator` | 广播所用的 MPI 通信域。 |
| `ierr` | `mpi_return_code` | MPI 广播返回码。 |

### `npol_abort_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `comm` | `communicator` | 要中止的 MPI 通信域。 |
| `message` | `message` | rank 0 输出的错误信息。 |
| `rank` | `rank` | 当前进程在指定通信域中的编号。 |
| `ierr` | `mpi_return_code` | MPI 查询或中止调用的返回码。 |

## [mpi_rmcdhf_csf_input.f90](../../src/lib/mpi/mpi_rmcdhf_csf_input.f90)

### `load_rmcdhf_descriptors_mpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `name` | `descriptor_path` | RMCDHF 的 CSF 描述符路径。 |
| `ncore` | `core_subshell_count` | 加载后返回的闭核子壳层数。 |
| `labels` | `block_labels` | 各 rank 已掌握的对称性块标签，用作一致性检查。 |
| `status` | `descriptor_status` | 读入、元数据一致性检查或数组分配的状态。 |
| `ierr` | `mpi_return_code` | MPI 广播返回码。 |
| `info` | `basis_counts` | rank 0 的 CSF 数、轨道数、电子数、块数和闭核子壳层数。 |
| `lengths` | `root_csf_counts_by_block` | rank 0 的各对称性块 CSF 数。 |
| `root_labels` | `root_block_labels` | rank 0 广播的对称性块标签。 |

## [mpi_isotope_input.f90](../../src/lib/mpi/mpi_isotope_input.f90)

### `setisompi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `isofile` | `isotope_file_path` | rank 0 读入的同位素文件路径。 |
| `j` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `i` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `k` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `nwin` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `ios` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `npy` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `naky` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `my` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `con` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `fkk` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `ey` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `pzy` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `dnorm` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `pa` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `qa` | — | 旧径向波函数加载器遗留声明，无引用。 |
| `ra` | — | 旧径向波函数加载器遗留声明，无引用。 |

## [mpi_radial_wavefunction_file.f90](../../src/lib/mpi/mpi_radial_wavefunction_file.f90)

### `setrwfmpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `name` | `wavefunction_path` | 待打开的 .rwf 文件路径。 |
| `ios` | `header_read_status` | 读取文件头时的 IOSTAT。 |
| `ierror` | `setup_status` | 依次接收 OPENFL 和 LODRWFMPI 的返回码。 |
| `g92rwf` | `file_header` | 文件中读出的六字符标记，应等于 G92RWF。 |

## [mpi_sparse_matrix_product.f90](../../src/lib/mpi/mpi_sparse_matrix_product.f90)

### `spicmvmpi`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `matrix_order` | 稀疏哈密顿量阶数，也是一列向量的长度。 |
| `m` | `vector_count` | 本次同时相乘的向量列数。 |
| `b` | `input_vectors` | 稀疏哈密顿量作用的输入向量列。 |
| `c` | `product_vectors` | 先累计各 rank 的局部贡献，再接收全局矩阵乘结果。 |
| `icol` | `column_index` | 当前 rank 负责的矩阵列索引。 |
| `ibeg` | `segment_start` | 当前列在本地稀疏数组中的起始位置。 |
| `iend` | `segment_end` | 本地稀疏段的末尾位置，指向对角元。 |
| `nelc` | `segment_element_count` | 当前存储段的元素数，含对角元。 |
| `iv` | `vector_index` | 输入和输出向量的列索引。 |
| `diag` | `partial_product` | 此前累计值加本列对角贡献的暂存值。 |
| `dl` | `off_diagonal_product` | DMERGE 返回的本列非对角乘积贡献。 |

## [mpi_change_directory.f90](../../src/lib/mpi/mpi_change_directory.f90)

### `sys_chdir`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `dir` | `directory_path` | 目标工作目录。 |
| `lendir` | `path_length` | 传给系统例程的路径字符数。 |
| `ierr` | `directory_status` | 切换目录操作的返回码。 |

## [mpi_current_directory.f90](../../src/lib/mpi/mpi_current_directory.f90)

### `sys_getwd`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `dir` | `directory_path` | 返回的当前工作目录。 |
| `lcwd` | `path_length` | 返回的工作目录有效字符数。 |
| `status` | `getcwd_status` | 查询工作目录的返回码。 |

## [mpi_make_directory.f90](../../src/lib/mpi/mpi_make_directory.f90)

### `sys_mkdir`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `dir` | `directory_path` | 要创建的目录路径。 |
| `lendir` | `path_length` | 传给系统例程的路径字符数。 |
| `ierr` | `directory_status` | 创建目录操作的返回码。 |
| `mkdir_error_message` | `mkdir_error_message` | 保留现有的创建目录失败消息缓冲。 |

## 外部变量与常量的本地名称

以下名称通过 `USE ... ONLY` 在用到它的源码内提供。定义模块和其他调用方继续使用原符号。
已明确表达含义的共享变量与容量常量沿用其定义名，例如 `orbital_count`、`hamiltonian_values` 和
`max_radial_points`。

| 使用处 | 定义模块 | 原名 | 本地名称 | 用途 |
| --- | --- | --- | --- | --- |
| `cpath` | `mpi_state` | `mpi_error_code` | `directory_status` | 传递目录操作的状态，与通信返回码区分。 |
| `load_csf_basis_mpi` | `blk_c` | `nblock` | `symmetry_block_count` | 对称性块数。 |
| `load_csf_basis_mpi` | `blk_c` | `ncfblk` | `cumulative_csf_counts` | CSF 按对称性块的累计计数。 |
| `load_rmcdhf_descriptors_mpi` | `syma_c` | `jpgg` | `packed_block_symmetries` | 打包的对称性块量子数编码。 |
| `load_rmcdhf_descriptors_mpi` | `syma_c` | `nblk0` | `max_symmetry_blocks` | 对称性块数组容量。 |
| `lodrwfmpi` | `debug_c` | `ldbpr` | `rmcdhf_debug_flags` | 第 3 项控制径向波函数加载诊断输出。 |
| `npol_init_mpi` | `npol_c` | `zz_coul` | `coulomb_potential_baseline` | 加上极化修正前的库仑势，采用 -rV 表示。 |
| `setisompi` | `def_c` | `emn` | `nuclear_mass_in_electron_masses` | 以电子质量为单位的核质量。 |
| `setisompi` | `nsmdat_c` | `sqn` | `nuclear_spin` | 核自旋。 |
| `setisompi` | `nsmdat_c` | `dmomnm` | `nuclear_magnetic_dipole_moment` | 核磁偶极矩。 |
| `setisompi` | `nsmdat_c` | `qmomb` | `nuclear_electric_quadrupole_moment` | 核电四极矩。 |

清理了 `lodrwfmpi` 中未使用的调试开关和部分网格、核模型、精度导入，`setisompi` 中未使用的
精度与正交化导入，以及 `setrwfmpi` 的输入和输出逻辑单元、`spicmvmpi` 的矩阵元数导入。

## 验证

- 按作用域核对 294 个原声明和迁移后的名称。按对应表替换标识符后，可执行语句与迁移前一致。
- GNU Fortran 15.3.0 与 Open MPI 的 Release 和 Debug 均构建通过，构建包含串行和 MPI 程序。
  109 项 Fortran 与集成测试中，Release 通过 109 项，Debug 通过 108 项。
- Debug 的 `rci_mpi_breit` 在迁移前后均于同一径向积分语句触发浮点异常。该测试在 Release 通过。
  [zkf.f90 第 102 行](../../src/apps/rci_mpi/evaluate_breit_hartree_z.f90#L102)
