# Davidson 本征求解库命名迁移

## 库、源码和模块

| 项目 | 原名 | 现名 |
| --- | --- | --- |
| 源码目录 | `src/lib/libdvd90/` | [`src/lib/davidson_eigensolver/`](../../src/lib/davidson_eigensolver/) |
| CMake 目标 | `dvd90` | `davidson_eigensolver` |
| 静态库 | `libdvd90.a` | `libdavidson_eigensolver.a` |
| 链接器参数 | `-ldvd90` | `-ldavidson_eigensolver` |
| 构建产物 | `<build>/lib/libdvd90.a` | `<build>/lib/libdavidson_eigensolver.a` |
| 安装产物 | `<prefix>/lib/libdvd90.a` | `<prefix>/lib/libdavidson_eigensolver.a` |
| 构建模块目录 | `<build>/modules/dvd90/` | `<build>/modules/davidson_eigensolver/` |
| 安装模块目录 | `<prefix>/lib/dvd90/` | `<prefix>/lib/davidson_eigensolver/` |

每个新文件的文件名主体和 Fortran 模块名相同。生成的 `.mod` 文件也使用该名称。

| 原文件 | 原模块 | 新模块 | 职责 |
| --- | --- | --- | --- |
| `adds.f90` | `adds_i` | [`davidson_projected_matrix`](../../src/lib/davidson_eigensolver/davidson_projected_matrix.f90) | 追加打包投影矩阵的列 |
| `dvdrvr.f90` | `dvdrvr_i` | [`davidson_iteration`](../../src/lib/davidson_eigensolver/davidson_iteration.f90) | Davidson 主迭代和重启控制 |
| `dvdson.f90` | `dvdson_i` | [`davidson_solver_control`](../../src/lib/davidson_eigensolver/davidson_solver_control.f90) | 目标选择、工作区和反向通信调度 |
| `gdvd.f90` | `gdvd_i` | [`davidson_eigensolver`](../../src/lib/davidson_eigensolver/davidson_eigensolver.f90) | 求解入口、矩阵–向量回调和对角预条件 |
| `iniest.f90` | `iniest_i` | [`davidson_sparse_initial_estimates`](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90) | 从稀疏对称性块生成初猜 |
| `initdvd.f90` | `initdvd_i` | [`davidson_initial_subspace`](../../src/lib/davidson_eigensolver/davidson_initial_subspace.f90) | 构造初始子空间 |
| `mgs_nrm.f90` | `mgs_nrm_i` | [`davidson_orthonormalization`](../../src/lib/davidson_eigensolver/davidson_orthonormalization.f90) | 改进 Gram–Schmidt 正交化与归一化 |
| `multbc.f90` | `multbc_i` | [`davidson_basis_transform`](../../src/lib/davidson_eigensolver/davidson_basis_transform.f90) | 原地变换基向量 |
| `newvec.f90` | `newvec_i` | [`davidson_residual_expansion`](../../src/lib/davidson_eigensolver/davidson_residual_expansion.f90) | 计算 Ritz 残差和扩展方向 |
| `ovflow.f90` | `ovflow_i` | [`davidson_subspace_restart`](../../src/lib/davidson_eigensolver/davidson_subspace_restart.f90) | 压缩子空间并重启 |
| `tstsel.f90` | `tstsel_i` | [`davidson_convergence`](../../src/lib/davidson_eigensolver/davidson_convergence.f90) | 检测收敛并排序待处理本征对 |

过程名和参数仍按现有接口使用，例如 `use davidson_eigensolver, only: gdvd, dvd_operator`。
外部 Fortran 调用方需要更新 `USE` 和链接名称，使用新的 `.mod` 文件重新编译和链接。
模块名是 Fortran 链接符号的一部分。

### 库和模块命名的验证

- 用 GNU Fortran 15.3.0 构建了串行及 MPI Release 目标，110/110 项 CTest 通过。
- 在独立的安装前缀核对了 `libdavidson_eigensolver.a` 和全部 11 份新 `.mod`。
  使用安装的模块和 `-ldavidson_eigensolver` 重新编译并链接 Davidson 测试程序后，
  两组夹具的数值结果和迭代计数与改名前逐位一致。

## 变量和形参

变量迁移以 2026-11-20 的 `2d4c721` 为基准，当时源码目录为 `src/lib/libdvd90`。下表按过程作用域列出
形参、局部变量及函数结果，每行对应一个作用域实体。表中的「原名」按 Fortran 不区分大小写的规则记录。

### 命名与接口约定

- `matrix_order` 是原矩阵阶数，`subspace_size` 是当前基向量数，`subspace_limit` 是其容量。
  `selected_pair_count` 只计目标态，而 `extreme_pair_count` 计从谱端到最远目标态必须覆盖的本征对。
  例如，从低谱端请求第 1 和第 3 态时，这两个计数分别为 2 和 3。
- `ritz_values` 和 `ritz_coefficients` 表示投影问题的本征值和本征向量。
  `operator_basis` 保存 `A B` 并复用尾部暂存空间，`projected_matrix` 保存打包后的 `Bᵀ A B`。
  `*_start` 表示 Fortran 一基起始位置，`*_offset` 表示相对偏移。
- 七个 `request_slots` 依次为当前批次向量数、输入位置、输出位置、Ritz 值表起点、整数工作区中
  待处理本征对列表起点、残差范数表起点、预条件位移表起点。后四项仅在预条件请求中设置。
  `initdvd` 和 `dvdrvr` 在第 2、3 槽写基向量列号，`dvdson` 将它们转换为 `real_workspace` 的首元素位置。
- 变量迁移阶段保留了过程和外部 BLAS/LAPACK 名称。过程参数的顺序、类型、数组布局和传递方式
  保持原有约定。仓库内相关调用使用位置参数。外部 Fortran 调用方若使用关键字参数，需要按表更新关键字。
- 在这 11 份 Fortran 源码中，`real64` 是从 `iso_fortran_env` 导入的 kind 常量，沿用标准名。
  库内没有本地声明的 `PARAMETER` 或 `DATA` 常量。数值字面量、协议编号和诊断字符串沿用原值。
  此目录的 `CMakeLists.txt` 定义静态库目标，没有本地 CMake 变量或常量声明。

本次登记 **277** 个作用域实体，重命名 **276** 个，核定沿用 **1** 个。

### [`davidson_projected_matrix.f90`](../../src/lib/davidson_eigensolver/davidson_projected_matrix.f90)

追加投影矩阵列。

#### `adds`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `lim` | `subspace_limit` | Davidson 子空间允许容纳的最大基向量数。 |
| `kpass` | `subspace_size` | 当前子空间已经积累的基向量数。 |
| `basis` | `basis_vectors` | 展平的全空间基向量矩阵，初始化和回写时也存放近似本征向量。 |
| `ab` | `operator_basis` | 展平的矩阵作用结果，即 A 乘基向量矩阵。尾部可暂存残差。 |
| `s` | `projected_matrix` | 按上三角打包的投影矩阵，等于 BᵀAB。 |
| `nncv` | `new_vector_count` | 本次需要补算投影列的新增基向量数。 |
| `idstart` | `operator_column_start` | 当前新增 A·B 列的起始位置。 |
| `isstart` | `packed_column_offset` | 当前投影列之前已经占用的打包元素数。 |
| `iv` | `new_vector_index` | 新增列在本次追加批次中的序号。 |
| `ibstart` | `basis_column_start` | 参与点积的基向量列起始位置。 |
| `ibv` | `basis_vector_index` | 参与点积的基向量序号，也是当前投影列的行号。 |
| `ss` | `projected_element` | 一次基向量点积得到的投影矩阵元。 |

### [`davidson_iteration.f90`](../../src/lib/davidson_eigensolver/davidson_iteration.f90)

Davidson 迭代、收敛与重启。

#### `dvdrvr`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `irc` | `operation_code` | 反向通信请求类型，0 就绪或结束，1 预条件，2 初始化矩阵乘，3 迭代矩阵乘。 |
| `irev` | `request_slots` | 反向通信的向量数及索引槽位。驱动层先写列号，入口层再转为工作数组位置。 |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `hiend` | `solve_high_end` | 是否从高能谱端求解。 |
| `lim` | `subspace_limit` | Davidson 子空间允许容纳的最大基向量数。 |
| `mblock` | `expansion_block_size` | 每轮允许加入的修正向量数上限。 |
| `noc` | `constraint_count` | 需要保持正交的约束向量数。 |
| `nume` | `extreme_pair_count` | 从选中的谱端到最远目标态所覆盖的本征对数。 |
| `niv` | `initial_vector_count` | 初始近似向量数。初始化过程可在补齐向量后更新它。 |
| `neig` | `selected_pair_count` | 实际请求的本征对数，稀疏选择时可以小于 extreme_pair_count。 |
| `iselec` | `selected_pair_indices` | 用户指定或由范围展开得到的本征对序号。 |
| `crite` | `eigenvalue_tolerance` | 相邻两轮本征值变化量的收敛阈值。 |
| `critc` | `coefficient_tolerance` | 新增基向量对应的 Ritz 系数的收敛阈值。 |
| `critr` | `residual_tolerance` | 本征对残差范数的收敛阈值。 |
| `maxiter` | `iteration_limit` | Davidson 迭代次数上限。 |
| `eigval` | `ritz_values` | 当前迭代中的 Ritz 值，即投影矩阵本征值。 |
| `basis` | `basis_vectors` | 展平的全空间基向量矩阵，初始化和回写时也存放近似本征向量。 |
| `orthobasis` | `basis_with_constraints` | 约束向量与当前基向量连续存放的正交化视图。 |
| `ab` | `operator_basis` | 展平的矩阵作用结果，即 A 乘基向量矩阵。尾部可暂存残差。 |
| `s` | `projected_matrix` | 按上三角打包的投影矩阵，等于 BᵀAB。 |
| `svec` | `ritz_coefficients` | 投影矩阵的本征向量，给出全空间 Ritz 向量在当前基下的系数。 |
| `incv` | `pending_pair_indices` | 待处理本征对序号。在 LAPACK 调用期间复用为失败向量序号缓冲。 |
| `icv` | `convergence_flags` | 按本征对序号记录收敛状态的整型标志。 |
| `oldval` | `previous_ritz_values` | 上轮 Ritz 值。`dvdrvr` 在返回前将其覆盖为绝对变化量。 |
| `nloops` | `iteration_count` | Davidson 迭代计数。 |
| `ierr` | `solver_status` | 求解器返回状态，含输入检查或迭代错误。 |
| `done` | `targets_converged` | 当前收敛判据是否允许结束迭代。 |
| `temps` | `lapack_packed_copy` | 投影矩阵的副本，供会覆盖输入的 LAPACK 求解使用。 |
| `scra1` | `real_scratch` | 复用的实数暂存区，用于 LAPACK、系数排序、残差范数、预条件位移和正交化。 |
| `iscra2` | `integer_scratch` | 复用的整数暂存区，用于 LAPACK 和候选本征对排序。 |
| `rest` | `restart_size` | 达到子空间上限时重启保留的向量数。 |
| `i` | `pair_position` | 遍历本征对或待处理索引列表时的序号。 |
| `kpass` | `subspace_size` | 当前子空间已有的基向量数。通过 SAVE 跨回调保留。 |
| `nncv` | `expansion_count` | 当前或上一轮扩展向量数。在筛选和残差计算之间暂存待处理本征对数。 |
| `ifind` | `requested_ritz_count` | 本轮要求 LAPACK 返回的投影本征对数。 |
| `nfound` | `computed_ritz_count` | LAPACK 报告实际求出的本征对数。 |
| `info` | `lapack_status` | LAPACK 的 INFO 返回值。 |
| `newstart` | `expansion_start` | 新增向量在展平数组中的起始位置。 |
| `tol` | `lapack_abs_tolerance` | 传给 LAPACK 的绝对容差，取 SLAMCH 的安全最小值。 |
| `first` | `first_iteration` | 历史首轮标志。现有实现赋值并保存，迭代中无读取。 |

### [`davidson_solver_control.f90`](../../src/lib/davidson_eigensolver/davidson_solver_control.f90)

目标选择、工作区布局与反向通信。

#### `dvdson`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `irc` | `operation_code` | 反向通信请求类型，0 就绪或结束，1 预条件，2 初始化矩阵乘，3 迭代矩阵乘。 |
| `irev` | `request_slots` | 反向通信的向量数及索引槽位。驱动层先写列号，入口层再转为工作数组位置。 |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `lim` | `subspace_limit` | Davidson 子空间允许容纳的最大基向量数。 |
| `noc` | `constraint_count` | 需要保持正交的约束向量数。 |
| `ilow` | `lowest_pair_index` | 请求的最低本征对序号，从 1 开始。 |
| `ihigh` | `highest_pair_index` | 请求的最高本征对序号，从 1 开始。 |
| `iselec` | `selected_pair_indices` | 用户指定或由范围展开得到的本征对序号。 |
| `niv` | `initial_vector_count` | 初始近似向量数。初始化过程可在补齐向量后更新它。 |
| `mblock` | `expansion_block_size` | 每轮允许加入的修正向量数上限。 |
| `maxiter` | `iteration_limit` | Davidson 迭代次数上限。 |
| `iwrsz` | `real_workspace_size` | 实数工作数组的元素数。 |
| `iiwsz` | `integer_workspace_size` | 整数工作数组的元素数。 |
| `nloops` | `iteration_count` | Davidson 迭代计数。 |
| `ierr` | `solver_status` | 求解器返回状态，含输入检查或迭代错误。 |
| `crite` | `eigenvalue_tolerance` | 相邻两轮本征值变化量的收敛阈值。 |
| `critc` | `coefficient_tolerance` | 新增基向量对应的 Ritz 系数的收敛阈值。 |
| `critr` | `residual_tolerance` | 本征对残差范数的收敛阈值。 |
| `hiend` | `solve_high_end` | 是否从高能谱端求解。 |
| `iwork` | `integer_workspace` | 整数工作数组，承载 LAPACK 暂存区及迭代索引和标志。 |
| `work` | `real_workspace` | 实数工作数组，承载初始基、迭代数据和输出。 |
| `neig` | `selected_pair_count` | 实际请求的本征对数，稀疏选择时可以小于 extreme_pair_count。 |
| `nume` | `extreme_pair_count` | 从选中的谱端到最远目标态所覆盖的本征对数。 |
| `ddot` | `unused_dot_product` | 未使用的实数历史声明。其他作用域导入的 BLAS DDOT 是独立的过程符号。 |
| `i` | `selection_position` | 选择列表位置，也暂存高谱端索引转换后的下一位置。 |
| `iortho` | `constrained_basis_start` | 实数工作区中含约束向量的基矩阵起始位置。 |
| `ibasis` | `basis_start` | 实数工作区中无约束部分的基矩阵起始位置。 |
| `ieigval` | `ritz_values_start` | 实数工作区中 Ritz 值的起始位置。 |
| `iab` | `operator_basis_start` | 实数工作区中 A 乘基矩阵的起始位置。 |
| `is` | `projected_matrix_start` | 实数工作区中打包投影矩阵的起始位置。 |
| `itemps` | `lapack_packed_copy_start` | 实数工作区中投影矩阵副本的起始位置。 |
| `isvec` | `ritz_coefficients_start` | 实数工作区中 Ritz 系数矩阵的起始位置。 |
| `iscra1` | `real_scratch_start` | 实数暂存区起始位置，预条件请求时起始位置存放残差范数。 |
| `ioldval` | `previous_ritz_values_start` | 上轮 Ritz 值暂存区的起始位置。 |
| `iscra2` | `integer_scratch_start` | 整数暂存区的起始位置。 |
| `iscra3` | `pending_indices_start` | 待处理本征对序号表在整数工作区中的起始位置。 |
| `iicv` | `convergence_flags_start` | 收敛标志在整数工作区中的起始位置。 |

### [`davidson_eigensolver.f90`](../../src/lib/davidson_eigensolver/davidson_eigensolver.f90)

对角预条件与矩阵–向量乘法回调。

#### `dvd_operator`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `m` | `vector_count` | 本次回调处理的列向量数。 |
| `b` | `input_vectors` | 输入列向量集合。 |
| `c` | `output_vectors` | 回调写入的矩阵–向量乘法结果。 |

#### `gdvd`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `lim` | `subspace_limit` | Davidson 子空间允许容纳的最大基向量数。 |
| `ilow` | `lowest_pair_index` | 请求的最低本征对序号，从 1 开始。 |
| `ihigh` | `highest_pair_index` | 请求的最高本征对序号，从 1 开始。 |
| `iselec` | `selected_pair_indices` | 用户指定或由范围展开得到的本征对序号。 |
| `niv` | `initial_vector_count` | 初始近似向量数。初始化过程可在补齐向量后更新它。 |
| `mblock` | `expansion_block_size` | 每轮允许加入的修正向量数上限。 |
| `maxiter` | `iteration_limit` | Davidson 迭代次数上限。 |
| `iwrsz` | `real_workspace_size` | 实数工作数组的元素数。 |
| `iiwsz` | `integer_workspace_size` | 整数工作数组的元素数。 |
| `nloops` | `iteration_count` | Davidson 迭代计数。 |
| `ierr` | `solver_status` | 求解器返回状态，含输入检查或迭代错误。 |
| `crite` | `eigenvalue_tolerance` | 相邻两轮本征值变化量的收敛阈值。 |
| `critc` | `coefficient_tolerance` | 新增基向量对应的 Ritz 系数的收敛阈值。 |
| `critr` | `residual_tolerance` | 本征对残差范数的收敛阈值。 |
| `hiend` | `solve_high_end` | 是否从高能谱端求解。 |
| `iwork` | `integer_workspace` | 整数工作数组，承载 LAPACK 暂存区及迭代索引和标志。 |
| `work` | `real_workspace` | 实数工作数组，承载初始基、迭代数据和输出。 |
| `noc` | `constraint_count` | 需要保持正交的约束向量数。 |
| `irc` | `operation_code` | 反向通信请求类型，0 就绪或结束，1 预条件，2 初始化矩阵乘，3 迭代矩阵乘。 |
| `matrix_vector_callback` | `matrix_vector_callback` | 执行矩阵–向量乘法的过程形参，沿用已有清晰名称。 |
| `diag` | `matrix_diagonal` | 矩阵对角元，执行每个预条件操作时临时减去相应的 Ritz 值。 |
| `ortho` | `orthogonality_tolerance` | 历史接口保留参数。此包装过程未读取或向下传递。 |
| `nmv` | `matrix_vector_product_count` | 计入各批次数量的矩阵–向量乘法总数。 |
| `irev` | `request_slots` | 求解器返回的操作向量数及工作数组索引。七个槽位的解释见上文。 |
| `nb` | `vector_count` | 当前反向通信请求中涉及的列向量数。 |
| `iw1` | `input_start` | 当前操作的输入在实数工作区中的首元素位置。 |
| `iw2` | `output_start` | 当前操作的输出在实数工作区中的首元素位置。 |
| `iw3` | `ritz_values_start` | 当前 Ritz 值表在实数工作区中的起始位置。 |
| `iiw` | `pending_indices_start` | 待处理本征对列表在整数工作区中的起始位置。 |
| `iin` | `residual_norms_start` | 残差范数表在实数工作区中的起始位置。 |
| `iw4` | `preconditioner_shifts_start` | 预条件位移表在实数工作区中的起始位置。 |
| `icur` | `vector_offset` | 当前操作向量之前的元素数，用于定位输入和输出列。 |
| `i` | `vector_index` | 批次中的向量序号。 |
| `indx` | `pair_offset` | 本征对序号减 1，用于在各个工作区切片中定位。 |
| `j` | `component_index` | 当前向量的分量序号，也是对角元下标。 |
| `value` | `preconditioning_value` | 此次对角预条件使用的 Ritz 值，求解器可在请求发出前调整该值的符号或位移。 |
| `ovalue` | `stored_shift` | 从第 7 槽指向的位移表读取的值，目前仅赋值，未参与公式。 |
| `rnorm` | `residual_norm` | 当前待处理本征对的残差范数，目前仅赋值，未参与公式。 |
| `epsil` | `preconditioner_shift` | 从残差范数表起点加 subspace_limit 读到的同一位移，目前仅赋值。 |

### [`davidson_sparse_initial_estimates.f90`](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90)

由对称性块生成初猜及分配失败诊断。

#### `iniest`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `niv` | `initial_vector_count` | 初始近似向量数。初始化过程可在补齐向量后更新它。 |
| `n` | `full_matrix_order` | 完整稀疏矩阵的阶数，也是输出初猜向量的长度。 |
| `nb` | `block_order` | 所选对称性块的阶数。 |
| `hmx` | `stored_matrix_values` | 稀疏矩阵的存储值。 |
| `jcol` | `column_end_offsets` | 各稀疏列末元素的位置，下标 0 记录起点边界。 |
| `irow` | `stored_row_indices` | 稀疏矩阵各已存元素的全矩阵行号。 |
| `basis` | `initial_estimates` | 输出工作区，依次存全空间初猜向量和对应本征值。 |
| `iblock` | `csf_block_indices` | 各全矩阵 CSF 所属的对称性块。 |
| `jblock` | `target_block` | 要提取并求初猜的对称性块编号。 |
| `ns` | `estimate_matrix_order` | 用于初猜的小矩阵阶数，当前上限 1000。 |
| `icount` | `block_column` | 已提取的块内列数，也是当前压缩列号。 |
| `j` | `column_index` | 当前全矩阵列号，回填时用于遍历初猜向量列。 |
| `istart` | `storage_index` | 提取时为稀疏元素位置，回填时为输出向量之前的元素数。 |
| `ish` | `row_index_shift` | 将全矩阵行号转换成所选块内行号时扣除的偏移。 |
| `ii` | `unused_index` | 未引用的历史整数声明，无法从现存代码给出更具体的用途。 |
| `nfound` | `computed_pair_count` | LAPACK 报告求出的本征对数。 |
| `info` | `lapack_status` | LAPACK 的 INFO 返回值。 |
| `ierr` | `status_code` | 分配状态或由 LAPACK 返回值换算的状态。 |
| `ic` | `block_row` | 把局部本征向量嵌回全空间时的块内行号。 |
| `iqiwork` | `lapack_integer_work_alias` | 指向 LAPACK 整数暂存区的历史指针别名。 |
| `iqif` | `failed_vector_indices_alias` | 指向失败本征向量序号表的历史指针别名。 |
| `iqjsh` | `block_to_full_index_alias` | 指向块内到全矩阵索引表的历史指针别名。 |
| `iwork` | `lapack_integer_work` | LAPACK 的整数暂存数组。 |
| `ifail` | `failed_vector_indices` | LAPACK 记录未成功求得的本征向量序号的数组。 |
| `jsh` | `block_to_full_index` | 选中块中各保留行列对应的全矩阵序号。 |
| `iqap` | `packed_estimate_matrix_alias` | 指向打包初猜矩阵的历史指针别名。 |
| `iqeig` | `initial_eigenvalues_alias` | 指向初猜本征值的历史指针别名。 |
| `iqvec` | `initial_block_vectors_alias` | 指向块内初猜本征向量的历史指针别名。 |
| `iqwork` | `lapack_real_work_alias` | 指向 LAPACK 实数暂存区的历史指针别名。 |
| `ap` | `packed_estimate_matrix` | 所选块的初猜矩阵，按下三角打包。 |
| `eigval` | `initial_eigenvalues` | 初猜小矩阵求得的本征值。 |
| `vec` | `initial_block_vectors` | 初猜小矩阵求得的本征向量，以块内坐标表示。 |
| `work` | `lapack_real_work` | LAPACK 的实数暂存数组。 |

#### `mem_fail`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `unit` | `output_unit` | 分配失败诊断的 Fortran 输出单元。 |
| `count` | `element_count` | 分配失败的数组元素数。 |
| `variable` | `variable_name` | 诊断中打印的数组标签。 |
| `status` | `allocation_status` | ALLOCATE 返回的状态值。 |

### [`davidson_initial_subspace.f90`](../../src/lib/davidson_eigensolver/davidson_initial_subspace.f90)

初始化子空间和第一份投影矩阵。

#### `initdvd`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `irc` | `operation_code` | 反向通信请求类型，0 就绪或结束，1 预条件，2 初始化矩阵乘，3 迭代矩阵乘。 |
| `irev` | `request_slots` | 反向通信的向量数及索引槽位。驱动层先写列号，入口层再转为工作数组位置。 |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `noc` | `constraint_count` | 需要保持正交的约束向量数。 |
| `niv` | `initial_vector_count` | 初始近似向量数。初始化过程可在补齐向量后更新它。 |
| `lim` | `subspace_limit` | Davidson 子空间允许容纳的最大基向量数。 |
| `hiend` | `solve_high_end` | 是否从高能谱端求解。 |
| `orthobasis` | `basis_with_constraints` | 约束向量与当前基向量连续存放的正交化视图。 |
| `basis` | `basis_vectors` | 展平的全空间基向量矩阵，初始化和回写时也存放近似本征向量。 |
| `ab` | `operator_basis` | 展平的矩阵作用结果，即 A 乘基向量矩阵。尾部可暂存残差。 |
| `s` | `projected_matrix` | 按上三角打包的投影矩阵，等于 BᵀAB。 |
| `kpass` | `subspace_size` | 构造首个投影矩阵时传给 adds 的已有列数，当前调用从 0 开始。 |
| `nume` | `minimum_initial_basis_size` | 初始化后需要达到的最小基向量数，由调用方传入 extreme_pair_count + 1。 |
| `scra1` | `orthogonalization_work` | 初始扩充向量正交化用的实数暂存区。 |
| `ist` | `first_batch_vector` | 当前初始化批次中第一个基向量的序号。 |
| `iend` | `last_basis_vector` | 当前已生成的最后一个基向量的序号。 |
| `new` | `batch_vector_count` | 当前初始化批次的向量数。 |

### [`davidson_orthonormalization.f90`](../../src/lib/davidson_eigensolver/davidson_orthonormalization.f90)

两遍改进 Gram–Schmidt 正交化和归一化。

#### `mgs_nrm`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `vector_length` | 每列向量的长度。 |
| `kp` | `existing_vector_count` | 正交化前已经存在的基向量及约束向量数。 |
| `new` | `new_vector_count` | 需要正交并归一的新增向量数。 |
| `scra` | `projection_coefficients` | 点积得到的投影系数暂存区。 |
| `b` | `vectors` | 包含已有向量和新增向量的展平矩阵。 |
| `i` | `pass_index` | 两遍改进 Gram–Schmidt 正交化的遍次。 |
| `newstart` | `new_vectors_start` | 新增向量在展平矩阵中的起始位置。 |
| `kcur` | `reference_vector_start` | 当前作为投影参考的向量起始位置。 |
| `k` | `reference_vector_index` | 当前作为投影参考的列向量序号。 |
| `jcur` | `target_vector_start` | 当前要减去投影的向量起始位置。 |
| `j` | `target_vector_index` | 当前要减去投影的列向量序号。 |
| `mm` | `component_offset` | 分量在其向量内的零基偏移。 |
| `dnm` | `vector_norm` | 当前被归一化向量的范数，先以平方值暂存，再开方。 |

### [`davidson_basis_transform.f90`](../../src/lib/davidson_eigensolver/davidson_basis_transform.f90)

原地右乘系数矩阵。

#### `multbc`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `row_count` | 要变换的矩阵行数。 |
| `k` | `input_column_count` | 变换前的列数，也是系数矩阵的行数。 |
| `m` | `output_column_count` | 变换后保留的列数，也是系数矩阵的列数。 |
| `c` | `coefficients` | 右乘的系数矩阵。 |
| `temp` | `row_workspace` | 单行变换结果的暂存向量。 |
| `b` | `vectors` | 原地右乘系数矩阵的展平向量矩阵。 |
| `irow` | `row_index` | 当前被变换的行号。 |

### [`davidson_residual_expansion.f90`](../../src/lib/davidson_eigensolver/davidson_residual_expansion.f90)

残差、预条件位移与扩展向量。

#### `newvec`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `nume` | `extreme_pair_count` | 从选中的谱端到最远目标态所覆盖的本征对数。 |
| `lim` | `subspace_limit` | Davidson 子空间允许容纳的最大基向量数。 |
| `mblock` | `expansion_block_size` | 每轮允许加入的修正向量数上限。 |
| `kpass` | `subspace_size` | 当前子空间已经积累的基向量数。 |
| `critr` | `residual_tolerance` | 本征对残差范数的收敛阈值。 |
| `incv` | `pending_pair_indices` | 待处理本征对序号。计算残差后保留实际纳入扩展的序号。 |
| `svec` | `ritz_coefficients` | 投影矩阵的本征向量，给出全空间 Ritz 向量在当前基下的系数。 |
| `eigval` | `ritz_values` | 当前迭代中的 Ritz 值，即投影矩阵本征值。 |
| `oldval` | `previous_ritz_values` | 上轮 Ritz 值。`dvdrvr` 在返回前将其覆盖为绝对变化量。 |
| `ab` | `operator_basis` | 展平的矩阵作用结果，即 A 乘基向量矩阵。尾部可暂存残差。 |
| `basis` | `basis_vectors` | 展平的全空间基向量矩阵，初始化和回写时也存放近似本征向量。 |
| `icv` | `convergence_flags` | 按本征对序号记录收敛状态的整型标志。 |
| `done` | `targets_converged` | 当前收敛判据是否允许结束迭代。 |
| `nncv` | `candidate_count` | 输入待检查本征对数，正常扩展返回时写成实际加入的向量数。 |
| `scra1` | `residual_norms` | 按本征对序号存放已计算的残差范数。 |
| `epsil` | `preconditioner_shifts` | 用于后续预条件的 Ritz 值位移，当前实现对已处理索引写零。 |
| `newstart` | `expansion_start` | 新向量在展平基和乘积数组中的起始位置。 |
| `nadded` | `added_vector_count` | 已接受的扩展向量数。 |
| `icvc` | `converged_candidate_count` | 残差范数达到阈值的候选本征对数。 |
| `limadd` | `expanded_size_limit` | 同时受剩余子空间容量和单轮扩展上限约束的子空间大小。 |
| `icur` | `expansion_vector_start` | 当前试探向量及其残差的起始位置。 |
| `i` | `candidate_position` | 当前候选在待处理列表中的位置。 |
| `indx` | `pair_index` | 当前候选对应的本征对序号。 |
| `sqres` | `squared_residual_norm` | 残差向量的平方范数。 |
| `gaplow` | `lower_ritz_gap` | 当前 Ritz 值与前一个值的间隔。 |
| `gap` | `nearest_ritz_gap` | 用于估计残差位移的最近 Ritz 值间隔。 |
| `gapup` | `upper_ritz_gap` | 当前 Ritz 值与后一个值的间隔。 |
| `dl` | `ritz_value_change` | 当前候选 Ritz 值相对上轮的绝对变化量。 |
| `res` | `residual_norm` | 当前候选本征对的残差范数。 |
| `rhseps` | `residual_shift` | 对预条件右端残差施加的沿 Ritz 向量方向的修正系数。 |

### [`davidson_subspace_restart.f90`](../../src/lib/davidson_eigensolver/davidson_subspace_restart.f90)

子空间压缩重启。

#### `ovflow`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `n` | `matrix_order` | 原实对称矩阵的阶数，也是一列全空间向量的长度。 |
| `kpass` | `subspace_size` | 当前子空间已经积累的基向量数。 |
| `basis` | `basis_vectors` | 展平的全空间基向量矩阵，初始化和回写时也存放近似本征向量。 |
| `ab` | `operator_basis` | 展平的矩阵作用结果，即 A 乘基向量矩阵。尾部可暂存残差。 |
| `s` | `projected_matrix` | 按上三角打包的投影矩阵，等于 BᵀAB。 |
| `svec` | `ritz_coefficients` | 投影矩阵的本征向量，给出全空间 Ritz 向量在当前基下的系数。 |
| `eigval` | `ritz_values` | 当前迭代中的 Ritz 值，即投影矩阵本征值。 |
| `nume` | `retained_vector_count` | 重启时保留的 Ritz 向量数。 |
| `scra1` | `row_workspace` | 原地变换基向量及其矩阵乘积时的逐行暂存区。 |
| `ind` | `packed_column_offset` | 重启后投影矩阵当前列之前的打包元素数。 |
| `icur` | `coefficient_column_offset` | 重置 Ritz 系数矩阵时当前列前的元素数。 |
| `i` | `pair_index` | 所保留 Ritz 本征对的序号。 |

### [`davidson_convergence.f90`](../../src/lib/davidson_eigensolver/davidson_convergence.f90)

目标收敛检测及候选排序。

#### `tstsel`

| 原名 | 迁移后 | 用途 |
| --- | --- | --- |
| `kpass` | `subspace_size` | 当前子空间已经积累的基向量数。 |
| `nume` | `extreme_pair_count` | 从选中的谱端到最远目标态所覆盖的本征对数。 |
| `neig` | `selected_pair_count` | 实际请求的本征对数，稀疏选择时可以小于 extreme_pair_count。 |
| `iselec` | `selected_pair_indices` | 用户指定或由范围展开得到的本征对序号。 |
| `svec` | `ritz_coefficients` | 投影矩阵的本征向量，给出全空间 Ritz 向量在当前基下的系数。 |
| `eigval` | `ritz_values` | 当前迭代中的 Ritz 值，即投影矩阵本征值。 |
| `icv` | `convergence_flags` | 按本征对序号记录收敛状态的整型标志。 |
| `crite` | `eigenvalue_tolerance` | 相邻两轮本征值变化量的收敛阈值。 |
| `critc` | `coefficient_tolerance` | 新增基向量对应的 Ritz 系数的收敛阈值。 |
| `oldval` | `previous_ritz_values` | 上轮 Ritz 值。`dvdrvr` 在返回前将其覆盖为绝对变化量。 |
| `incv` | `pending_pair_indices` | 待处理本征对序号。计算残差后保留实际纳入扩展的序号。 |
| `tstsel` | `targets_converged` | 函数结果变量，使用显式 RESULT 名称，过程仍为 tstsel。 |
| `nncv` | `recent_or_pending_count` | 入口为最近新增的向量数，出口为系数尚未收敛的目标本征对数。 |
| `rowlast` | `recent_coefficient_maxima` | 各未收敛目标在最近新增基向量上的绝对系数最大值，用于排序。 |
| `ind` | `unsorted_pair_indices` | 排序前与系数最大值对应的本征对序号。 |
| `nnce` | `unconverged_value_count` | 相邻轮次的 Ritz 值变化未满足容差的目标态数。 |
| `i` | `target_position` | 目标列表或当前排序列表中的位置。 |
| `ival` | `pair_index` | 当前目标对应的本征对序号。 |
| `icnt` | `pending_count` | 最新基向量的系数未达到阈值的目标数。 |
| `icur` | `last_coefficient_index` | 当前目标在展平 Ritz 系数矩阵中的最后一个系数下标。 |
| `l` | `coefficient_offset` | 相对于当前目标最后一个系数向前回溯的距离。 |
| `indx` | `relative_max_index` | IDAMAX 返回的绝对最大系数在未排序尾部的相对位置。 |
| `itemp` | `saved_pair_index` | 交换排序位置时暂存的本征对序号。 |
| `tmax` | `max_recent_coefficient` | 当前目标最近新增基向量上的绝对系数最大值。 |
| `temp` | `saved_coefficient_maximum` | 交换排序位置时暂存的系数最大值。 |
| `done` | `coefficients_converged` | 所有目标是否已通过最新基向量系数的收敛检查。 |

### 变量迁移的验证

- 逐作用域比对了 277 个声明及函数结果与对应表。按作用域反向映射新标识符后，11 份源码的声明和执行语句
  的词元与基线一致。显式 `RESULT` 子句为原有函数结果声明了独立名称。扫描到的 36 处 `CALL` 均使用位置参数。
- 使用 GNU Fortran 15.3.0 构建串行及 MPI 目标。Release 的 CTest 通过 110/110，Check 通过 109/110，
  两个 Davidson 专项用例在两种配置中均通过。
- 对照改名前后的 Release 库，耦合三对角矩阵与提前收敛矩阵两组结果的初猜、本征向量、本征值、
  本征值变化量、残差范数及迭代和矩阵乘法计数逐位一致。对照使用仓库的 Davidson 夹具。
  [夹具源码](../../test/f2018/davidson_contracts.F90)
- Check 唯一失败项 `rci_mpi_breit` 在进入本征求解之前的 Breit 积分阶段触发 `SIGFPE`。
  使用同一 Check 构建产物重新链接改名前的 `libdvd90` 后，同一用例在相同位置复现。
  [浮点异常位置](../../src/apps/rci_mpi/evaluate_breit_hartree_z.f90#L102)
