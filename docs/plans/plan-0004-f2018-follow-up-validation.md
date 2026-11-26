# Suite F2018 现代化实施与验收计划

**状态：P1 分批实施中。** 2026-11-26 按 Suite `3.4dev1` 的 `f211a0a35a957fab85c65a8a514d779183d0d59e` 复核。P0 提交 `88cb04a` 完成 RCI 续算边界、单 rank 失败传播及 V3 首轮对照。P1 的 MPI 绑定提交为 `5d77d93`，LP64 配置提交 `f493c99` 是稀疏矩阵批次的改造前基准。后续批次已实施基组与活动系数所有权、混合流写端、Davidson 残差判定和 Breit 故障传播，并取得四种 RMCDHF 的小型多根基准对照。

范围覆盖 `src/apps` 内 25 个 Fortran 应用目标，以及 `src/lib` 内 7 个库目录，含 Rust 描述符静态库及其 Fortran 调用端。OpenMP 的现行应用入口是 `jj2lsj_omp`。普通 RMCDHF 保留四个变体，各自有版本基线。目标范围以 [CMake 入口](../../CMakeLists.txt#L158)、[jj2lsj 构建定义](../../src/apps/jj2lsj/CMakeLists.txt)和[退役范围](../implemented/adr-0002-retire-strict-and-diagnostic-rmcdhf.md)为准。

## 当前基线和优先级

| 已核对的事实 | 对实施的影响 |
| --- | --- |
| 现有代码已采用显式过程接口、标准 kind 和严格构建配置。生产 `f2018-release` 默认用 `-O2`，Check/Audit 开启 `-std=f2018 -pedantic-errors` 和隐式接口诊断。GNU 16.2.1 下已有 `0f35673` 的 24 配置通过记录。 | 在基准和候选提交各记录实际编译命令、依赖库和测试清单。标准符合性、运行时检查、警告变化及数值回归各有验收结果。[构建配置](../../cmake/FortranProfiles.cmake) · [验收矩阵](../../test/f2018/README.md) |
| 当前 CSF 输入是 V3 Parquet 描述符及 TOML 头文件。波函数和混合系数分别要求 `RAD--WFN` 和 `CSF--MIX`，双正交输出也使用当前流格式。RCI 的 `.res` 是另一套续算协议。 | 固定输入的顺序、块身份和摘要。比较跨格式历史数据时按物理字段解码，比较本次同格式、同运算顺序的改造时可加字节对照。[文件契约](../../README_ZH.md) · [RCI 续算写端](../../src/apps/rci/read_ci_options.f90#L211) |
| 必需测试已覆盖 BLAS 调用、Davidson 小矩阵、MPI 消息、部分缓存和错误路径。`state_workflow` 已要求在两个对称块选根，检查根编号、缺根失败及已有输出。 | 在现有测试上扩展真实多根轨迹、态匹配与性质对照。历史上的 Zn I 三轮、Cd I 一轮 MCP 对照用于核验指定迭代和 I/O 行为，正常收敛需要独立算例。[已注册清单](../../test/f2018/required.json) · [选根用例](../../test/test_state_workflow.py#L197) · [MCP 对照](plan-0005-rmcdhf-mcp-production-io.md) |

本文的 **P0** 表示已定位且影响正常生产入口的错误处理缺口，以及验证修复所需的最小基线。**P1** 表示按调用链推进的绑定、接口和科学数值验收。**P2** 表示有实际消费者、所有权冲突或重复维护成本时再实施的设计改进。历史 GNU 15 方案的 P0–P9 是当时的阶段编号。F2018 允许模块状态、显式 `SAVE`、`POINTER` 和 assumed-size 形参，是否重构要结合有效输入、别名、生存期、并发及可检验性判断。[兼容性实施记录](../implemented/design-0004-f2018-compatibility.md) · [历史方案](../archive/plan-0002-f2018-migration-gnu15.md)

## 第一批实施

### P0 · 在当前入口建立首轮对照

**原有缺口。** 独立首轮比较器原先只准备 `isodata` 和 `rcsf.inp`。当前初猜从 header 取轨道表，角系数和 SCF 从 V3 文件对取基组。更新后的比较器生成一次输入，两套构建使用相同字节，分别保存输入、输出及摘要。`f2018_first_sweep_driver` 覆盖比较器本身，实际二进制对照使用独立入口。[比较器](../../test/f2018/compare_first_sweep.py) · [CTest 注册](../../test/f2018/CMakeLists.txt#L83) · [V3 夹具生成器](../../test/native_fixture.py)

**实施。** 用现有夹具生成器在独立目录生成配对的描述符和 header，给两侧提供相同的输入字节、状态名和交互应答。更新文件清单，既识别生成的输入，也检查输出完整。用相同编译器、数学库、优化选项及线程数分别构建基准和候选，存下输入、可执行文件摘要、标准输出和错误输出。先用当前源码的两个独立构建验证比较器能稳定识别一致结果，再以坏签名、缺块或改动一个科学值验证它能报告差异。历史版需要旧格式时，单独固定参考输入转换和字段映射。

**交付与退出。** 对 `rmcdhf` 和 `rmcdhf_mem` 分别取得可重放的一轮结果。核对 MCP 块与记录、轨道身份及 P/Q、ASF 根、能量和系数。明确预期停因为迭代上限。较大的真实多根算例在涉及相应内核改动前冻结，首轮软件回归可以独立完成。[现有工作流夹具](../../test/test_state_workflow.py) · [首轮诊断](../../test/f2018/compare_first_sweep.py#L110)

### P0 · 修正 RCI 续算记录的读取边界

**修复前。** `rci` 与 `rci_mpi` 的 `read_ci_restart` 在一条顺序无格式 `READ` 中先取得 `radial_point_count`，随即以它为上界读入核势数组，之后才检查 `max_radial_points`。核模型参数个数和每轨道尾点也在验证前用于读入。GNU Fortran Check 已复现这些坏记录造成的运行时错误。当前两种程序通过公共读取器先验证计数、回退记录并用已验证的边界读取。矩阵续算扫描也会核验块头、行元素数及列索引。[公共读取器](../../src/lib/core/rci_restart_records.f90) · [核参数容量](../../src/lib/grasp_shared/nuclear_model_parameters_state.f90) · [径向容量](../../src/lib/grasp_shared/cal_constants.f90)

1. 从写端列出每条记录的字段、计数、目标数组和后继使用者。用实际写端生成的小型合法文件作为正例。点核的参数个数为 0，现有有限核模型为 2。核对活动网格、含 10 个填充点的工作容量、每轨道尾点，以及 CSF 和对称块计数。[模型输入](../../src/lib/core/isotope_data_records.f90) · [串行写端](../../src/apps/rci/read_ci_options.f90#L215)
2. 在依赖计数的读取前验证取值和目标形状。对“计数与数组同记录”的情况，可先仅读计数，检查成功后 `BACKSPACE` 回到该记录再完整读取，对探读、定位和重读分别检查 `iostat`。顺序无格式文件中两次普通 `READ` 会消费两条记录，须用既有合法文件验明记录位置。尽量将记录解析供串行和 MPI 共用，由入口决定退出方式。
3. 按文件名、记录及适用的块或轨道报告首个错误。用合法边界值、越界计数、负计数、截断记录和非法轨道尾点做回归。正例做同模型冷启动与续算的态、能量和系数对照。负例要求程序主动给出可定位的非零失败，避免仅依赖 Debug 的运行时检查。
4. 单列基组身份问题。现有头记录可核对电子、CSF、轨道和块数，同维度文件仍可能对应不同的 CSF 排列。先追踪完整写入、读取与分片信息，再评估是否增加摘要或版本。若需改动磁盘协议，另列读写端覆盖、既有文件处置及同维异基组回归，把该项与读取边界修复独立提交。

### P0 · 处理 RCI 的单 rank 分片失败

**修复前。** `rci_mpi` 通过 root 广播续算选择，各 rank 随后打开自己的 `rciNNN.res`。某一片缺失、打开失败或签名错误时，该 rank 可单独进入带 barrier 的 `finalize_mpi_runtime`。用 Open MPI 5.0.11 复现 2 ranks 下 root 或 worker 单片缺失、签名错误都达到 8 秒测试超时。当前分片解析失败由 `mpi_io_errors` 报告记录和文件后终止作业。各 rank 只有部分矩阵写完时，先汇总是否存在已写的块头，再进入相同的统计归约分支。[分片入口](../../src/apps/rci_mpi/connect_ci_restart_file.f90) · [MPI I/O 退出](../../src/lib/mpi/mpi_io_errors.f90) · [矩阵续算](../../src/apps/rci_mpi/assemble_ci_matrix.f90)

**实施。** 在所有 rank 均能到达的阶段汇总本地状态，或在不可恢复的本地 I/O 错误上复用已有 `mpi_io_errors` 终止通路。保证出错 rank 能报告文件、rank 和阶段。沿 RCI 的头记录与矩阵记录调用链排查其他依赖集体到达的错误出口。正常完成仍经应用的集体收尾路径。[MPI I/O 错误处理](../../src/lib/mpi/mpi_io_errors.f90) · [续算矩阵乘](../../src/apps/rci_mpi/multiply_restart_hamiltonian.f90)

**交付与退出。** `rci_restart_serial` 和 `rci_restart_mpi_{1,2,4}` 调用实际程序生成合法文件并续算。四 CSF、两根夹具比较冷启动、完整续算和从干净记录边界补算的结果。串行和 2-rank 用例另覆盖缺文件、坏签名、坏计数及短记录，MPI 分别只破坏 root 或 worker 的一片。测试要求诊断含分片及读取阶段、命令在超时前非零退出。[回归入口](../../test/test_rci_restart.py) · [运行说明](../../test/f2018/README.md#rci-restart-records)

本机 GNU Fortran 16.2.1、OpenBLAS 与 Open MPI 5.0.11 的首批验收记录如下。

- 两套独立 Check 构建在改动 RCI 前完成 V3 单轮基线。RCI 改动后复跑，`rmcdhf` 和 `rmcdhf_mem` 各自的七份二进制文件、科学文本和规范化诊断仍一致，停因为一轮迭代上限。
- 串行以及 MPI 1、2、4 ranks 的四 CSF、两根夹具中，完整续算和文件尾补算与各自冷启动的最大能量差、相位对齐后的最大系数差均为 0。修复前写出的串行及双 rank 文件可由修复后程序续算。串行另外测了核模型参数为空的点核记录。
- 串行及 MPI 双 rank 的坏记录在 Check 和 Release 均非零退出并报告位置。合法续算仍使用原有顺序无格式记录布局。显式设置续算 ACCY 时，网格脚本修改串行和 MPI 共用的赋值位置。[网格脚本说明](../../scripts/README.md)

## P1 MPI 绑定、计数及分块

Suite 的生产和测试 MPI 调用统一使用 `mpi_f08`。`mpi_state` 继续提供共享进程状态并转发实际 MPI 接口。CMake 在启用 MPI 时检查 `MPI_Fortran_HAVE_F08_MODULE`，通过 `MPI::MPI_Fortran` 向消费者传递编译和链接依赖。核极化入口使用 `type(MPI_Comm)`，连续字节类型使用 `type(MPI_Datatype)`，点对点消息使用 `MPI_Status` 或明确的 `MPI_STATUS_IGNORE`。rank、tag 和错误码使用默认整数。[MPI 构建入口](../../CMakeLists.txt#L170) · [运行时接口](../../src/lib/mpi/mpi_runtime.f90) · [MPI Fortran 绑定](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node466.htm)

| 调用链 | 计数与缓冲约束 |
| --- | --- |
| CSF 装载与广播 | 保留每批最多 8192 个 CSF 的上限，再按最宽列和默认 MPI 元素数限额缩小批次。数组形状、可选字段和批次大小先在各 rank 对齐。通过 Fortran 的 `contiguous` 形参传入连续批次，跨步实参由编译器打包并回写。零列不发送 payload。[描述符广播](../../src/lib/mpi/mpi_csf_collectives.f90) |
| 两个 MPI RMCDHF 变体 | 直接系数以各 rank 的实际元素数执行 `MPI_Allgatherv`，用宽整数检查累计元素位移后转成默认整数界面。空贡献发送零个元素。保留按 rank 合并的顺序与径向求解器既有的有效系数数目。[RMCDHF 汇总](../../src/apps/rmcdhf_mpi/improve_orbital_mpi.f90) · [位移校验](../../src/lib/mpi/mpi_element_counts.f90) |
| 稀疏向量、双正交、密度和超精细矩阵 | 归约前以 `MPI_COUNT_KIND` 取得数组长度或计算维度乘积，再通过公共范围检查转换。现有密集初值与 RCI 稠密分支的阶数分别受求解器门槛约束。短标量、轨道和块元数据沿用相应 MPI 数据类型。[稀疏矩阵乘](../../src/lib/mpi/mpi_sparse_matrix_product.f90) · [双正交归约](../../src/apps/rbiotransform_mpi/transform_ci_coefficients.f90) · [密度归约](../../src/apps/rdensity/density_parallel_mpi.f90) |
| 文件及数值 ABI | `MPI_COUNT_KIND` 用于内存中的计数计算。输入输出文件的索引、记录和偏移使用各自格式规定的 kind。BLAS/LAPACK 仍通过 Suite 的 LP64 显式接口调用，MPI datatype 句柄不参与这些数值调用。[MCP 目录](../../src/lib/mpi/mpi_indexed_mcp.f90) · [BLAS/LAPACK 接口](../../src/lib/grasp_shared/blas_lapack_interfaces.f90) |

GNU Fortran 16.2.1 与本机 Open MPI 5.0.11 提供 `mpi_f08`，FindMPI 报告 MPI 标准版本 3.1。试编译仅广播一个元素的 `1_MPI_COUNT_KIND` 形式时，编译器报告没有匹配的 `MPI_Bcast` 重载。本次验收覆盖默认整数 collective 和有界分块路径。MPI-4 大计数重载须用提供该接口的部署库验证。[大计数标准接口](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node499.htm)

### 验收范围

- 1、2、4 ranks 用类型化字节及宽整数消息、点对点状态、可变长度汇总和解析求和验证接口。编译测试接受 `MPI_Datatype`，拒绝传入整数和 `MPI_Comm`。[消息契约](../../test/f2018/mpi_contracts.f90) · [类型检查](../../test/f2018/check_interface_rejections.py)
- 描述符测试覆盖一列一批、跨步缓冲区、部分最后一批、零列、一个 rank 出错与形状不一致。实际 MCP 稀疏写入和读回在 2、4 ranks 下验证一份有列和其余空分片。[描述符用例](../../test/descriptor_collectives.f90) · [MCP 空分片用例](../../test/f2018/sparse_empty_rank.f90)
- 2、4 ranks 的窄化回归令 rank 1 提供超限计数，要求在调用 payload collective 前非零退出，并核对源诊断。实际 RCI 续算另覆盖 1、2、4 ranks 合法分片及双 rank 单分片错误。[计数失败检查](../../test/f2018/check_mpi_count_failure.py) · [RCI 回归](../../test/test_rci_restart.py)

本机 `f2018-check` 构建在 MPI 与串行目标均启用时通过 19 项 MPI 专项和另外 108 项 CTest，关闭 MPI 时另核对配置的目标范围。在相同 GNU Fortran、OpenBLAS、V3 输入及同 rank 布局下，普通和内存 MPI RMCDHF 各运行一轮。1、2、4 ranks 的 `.wfn` 和 `.mix` 与各自 P0 基线逐字节一致。RMCDHF 导入来源校验匹配 248 个文件。[来源快照](../reference/rmcdhf_mpi_f08_snapshot.json)

## P1 BLAS/LAPACK LP64 接口与配置验收

Suite 通过 `grasp_shared` 的显式 Fortran 接口使用 BLAS/LAPACK。RMCDHF 的矩阵求解调用 `DCOPY/DDOT/DSCAL`，Davidson 使用 `DGEMV` 和返回整数索引的 `IDAMAX`。RCI 和 Davidson 的压缩矩阵求解使用 `DSPEVX`，其矩阵阶数、根索引、前导维度、返回根数、`INFO` 和整数工作区均为默认整数。配置探针确认该 kind 为 32 位，`SLAMCH` 和 `DLAMCH` 分别返回 32 位与 64 位实数。MPI 大计数及文件偏移有独立的宽度约定。[显式接口](../../src/lib/grasp_shared/blas_lapack_interfaces.f90) · [Davidson 残差扩展](../../src/lib/davidson_eigensolver/davidson_residual_expansion.f90) · [RCI MPI 特征求解](../../src/apps/rci_mpi/diagonalize_hamiltonian_block.f90)

CMake 固定以 `BLA_SIZEOF_INTEGER=4` 搜索，并在发现库后运行探针。BLAS 探针单独链接 `BLAS::BLAS` 并执行有界的 `DGEMV` 和 `DDOT` 对照。LAPACK 探针和正式目标同用 `LAPACK::LAPACK` 的传递链接顺序，验证零阶 `DSPEVX` 只写入 `M` 和 `INFO` 各自的 32 位槽，再检查浮点返回与 BLAS 调用。LAPACK 目标已传递 BLAS 依赖。两项探针每次配置均重新编译并运行，替换库后会重新验收。它们的编译和运行记录保存在构建目录的 `CMakeFiles/grasp-linalg/`。[配置验收](../../cmake/BlasLapackLP64.cmake) · [运行探针](../../cmake/blas_lapack_lp64_probe.F90) · [CMake 库宽度选项](https://cmake.org/cmake/help/v3.31/module/FindBLAS.html)

| 验收项 | 本机 GNU Fortran 16.2.1 / CMake 4.5.1 结果 |
| --- | --- |
| 正常配置及运行链接 | Check 和 Release 都选中 LP64 OpenBLAS 并通过两项配置探针。RCI 和线性代数测试的实际动态加载检查解析到 LP64 的 `libopenblas.so.0`。运行回归覆盖指针子数组和步长、`DDOT` 返回值、`DSPEVX` 本征值、残差及正交性。[生产接口回归](../../test/f2018/linear_algebra_contracts.f90) |
| 请求 ILP64 或调用者整数为 64 位 | 显式给出 `BLA_SIZEOF_INTEGER=8` 时，CMake 在库搜索前报错。构造的 64 位默认整数调用者也在探针中被拒绝。[配置回归](../../test/f2018/check_blas_lapack_config.py) |
| 实际选中同名符号的 ILP64 OpenBLAS | 本机两种库均为 OpenBLAS 0.3.34，其配置接口确认测试的 64 位库启用了 `USE64BITINT`。单独替换所选 BLAS 时，`DGEMV` 探针拒绝该库。保留 LP64 BLAS、只替换所选 LAPACK 时，BLAS 探针通过，`DSPEVX` 的输出槽检查拒绝 ILP64 LAPACK。回归用可独立编译的同符号库覆盖两个路径，并复用已成功配置的目录验证重新检查。 |

本机 Check 构建通过 9 项专项和其余 119 项 CTest，涵盖 RCI 串行及 MPI 1、2、4 ranks 的续算和 Davidson 求解。禁用 MPI 的 Release 构建通过配置验收与线性代数两项回归。

交叉配置通过 `CMAKE_CROSSCOMPILING_EMULATOR` 执行探针。配置和部署应使用相同的库与加载环境。更换编译器、供应商库或部署环境后需重新配置并运行线性代数回归。

### P1 · 稀疏矩阵存储与数组所有权

RCI 的 `.res` 按列存储，每列元素数和 CSF 索引为默认整数。串行和 MPI 原先用默认整数累加内存列尾，即使读入数组按 `int64` 的总元素数分配。串行 RCI 还在进入小矩阵分支前计算了默认整数的三角长度。现在累计位置从 RCI 装载一直以 `int64` 传到矩阵乘和 Davidson 初始子空间提取。小矩阵的三角长度只在阶数检查之后计算。RMCDHF 在顺序 MCP 及 MPI 索引 MCP 读入后扩展列尾，内存版从各自的 MCP 缓存扩展。MCP 的 32 位位置字段、RCI 每列记录和 LP64 调用继续使用各自的整数约定。[共享存储](../../src/lib/grasp_shared/sparse_hamiltonian_storage.f90) · [串行 RCI 装载](../../src/apps/rci/diagonalize_hamiltonian_block.f90) · [MPI RCI 装载](../../src/apps/rci_mpi/diagonalize_hamiltonian_block.f90) · [MPI MCP 读入](../../src/apps/rmcdhf_mpi/solve_hamiltonian_blocks_mpi.f90)

矩阵值、列尾和 CSF 索引由模块中的 `allocatable` 持有。RCI 装配期间按列容量分配，求解时按所选存储方式分配。RMCDHF 按块分配，轨道势例程在其独立阶段使用索引数组。调用链上没有这三个数组的原地扩容或指针别名。普通串行及两个内存版的跳块分支原先漏释放矩阵值，现将值数组分配放到跳过判定之后。串行 RCI 的内存估算及 MPI RMCDHF 的工作区预留已计入每个列尾 8 字节。通用 `ralloc` 的显式下界整数和字节重载按新旧下标区间的交集复制，覆盖空交集。GNU 16.2.1 在旧实现中复现了从 `-2:1` 扩到 `-3:3` 时读取旧数组下标 `-3` 的越界。[数组回归](../../test/f2018/array_memory.f90) · [跳块及重用回归](../../test/f2018/ci_damping.F90) · [RMCDHF 来源清单](../reference/rmcdhf_sparse_storage_snapshot.json)

本地用 GNU Fortran 16.2.1、LP64 OpenBLAS 和 Open MPI 5.0.11 构建候选工作树，在 `f2018-check` 下验证。串行构建通过 83 项 CTest。MPI 构建首轮通过 126 项，更新来源清单后复跑两个来源检查通过，使用较短临时路径重跑双基组流程也通过，共 129 项。RCI 续算和生产块求解的分配、跳块及重用覆盖 MPI 1、2、4 ranks。另有稀疏矩阵乘与空分区用例。和 `f493c99` 的串行固定首轮对照中，两种求解器各七份二进制输出逐字节一致，科学文本一致。内存版的主机空闲内存日志随运行时变化。两个 MPI 变体用同一套双对称块输入分别在 1、2、4 ranks 运行，改造前后均由求解器在第 2 轮结束。各布局的能量输出、混合系数和径向波函数逐字节一致，普通 MPI 版的每 rank 矩阵及径向追踪也逐字节一致。构建目录和日志在工作区的 `temp/tests/suite-p1-sparse/`。本批未运行超过 32 位累计位置的完整大矩阵或全配置验收矩阵。

### P1 · 基组与活动系数、混合流及多根回归

本批在上一批稀疏矩阵工作树基础上验证。活动与左侧基组的占据和耦合数组使用 `allocatable`，完整加载和合并使用 `move_alloc` 转移所有权。活动本征值、混合系数和 RCI 预计算系数使用原生分配与释放。描述符的 Rust 句柄由 Fortran 类型内部管理，关闭后清除摘要、块目录及轨道表。用两套不同 V3 基组交替读取，核对一侧失败后另一侧的已知数据和摘要，覆盖重复关闭、重开、完整基组替换以及打开失败时保留原基组。[完整基组加载](../../src/lib/core/csf_full_basis_reader.f90) · [描述符回归](../../test/descriptor_lifecycle.f90) · [来源清单](../reference/rmcdhf_basis_ownership_snapshot.json)

RCI 用两条合法的顺序记录验证累计位置跨过 `2^31-1`，MPI MCP 用小分片分别检查 `2^31-1` 与 `2^31` 的拒绝路径。混合流按块顺序写入，每一轮重新建立不完整流。刷出并关闭载荷后，核对磁盘长度，再写完成签名。固定写入上限复现了短写仍产生完成签名的问题，修复后的首轮与重复写入均保留不完整签名。MPI 根进程出错时调用现有的 MPI 终止回调。[RCI 计数夹具](../../test/f2018/sparse_offset_counts.f90) · [MPI MCP 边界](../../test/test_rmcdhf_mcp_paths.py) · [混合流写端](../../src/lib/core/mixing_file_format.f90) · [混合流回归](../../test/f2018/check_mixing_output.py)

解析双根夹具的能隙为 `1e-7`，有限试探空间触发 Davidson 重启。改造前的两个残差范数约 `2.85e-8`，求解器已返回成功。能量和系数的变化现用于确定扩展优先级，根的接受由残差决定。修复后两个残差范数约 `6.6e-13`，通过正交性、近简并子空间和重用返回向量的检查。Breit 小夹具对照直接项、交换项及相位，非法轨道组合和零张量秩分支输出诊断并以非零状态退出，MPI 单个工作 rank 的失败触发全作业终止。[Davidson 回归](../../test/f2018/davidson_contracts.F90) · [Breit 回归](../../test/f2018/breit_radial_contracts.F90) · [RCI MPI 接线](../../src/apps/rci_mpi/rci_mpi_main.f90)

四种 RMCDHF 均与各自的上一批构建对照，使用同一套 Ni 核参数、双对称块和初始轨道，每个布局共用一份 MCP 输入。偶块选根 1、2，奇块选根 1，等权并仅更新 `1s`。固定首轮均因迭代上限结束，独立运行中两种串行程序均在第 3 轮结束，两种 2-rank MPI 程序均在第 2 轮结束。各自对照的能量轨迹、混合流和径向波函数一致。固定轨道 RCI 的小稠密块由实际 `.res` 重建矩阵，核对冷启动、完整续算和部分续算的两个根。串行及 MPI 1、2、4 ranks 的最大归一化残差为 `9.548e-17`，最大正交误差为 `3.331e-16`。[RCI 残差回归](../../test/test_rci_restart.py)

GNU Fortran 16.2.1、LP64 OpenBLAS 和 Open MPI 5.0.11 的 `f2018-check` 构建通过串行 87 项及 MPI 135 项 CTest。用 16 个轨道和 8,193 个 CSF 另行核对 Fortran 读取边界，检查第 8,192 和 8,193 行的占据与次序，以及末行耦合码，读取进程峰值 RSS 为 143,276 KiB。构建及日志位于工作区 `temp/tests/suite-p1-follow-up/`。代表性大基组的峰值内存、真实多根生产模型与全配置矩阵仍需验收。

## 后续批次与共享库

先用第一批夹具保护读写和失败语义，再迁移通信接口。改动共用矩阵、轨道或角系数之前，先取得相应消费者的固定输入数值基线。单个子项目可按自身依赖推进。以下路径以实现库为单位；描述符的 Rust 端和 Fortran 调用端共同验收。

| 位置与优先级 | 实施范围和已有基础 | 验收证据 |
| --- | --- | --- |
| `mpi` · P1 | 本机迁移覆盖 communicator、datatype、status、常规元素计数和已有分块调用。部署库升级时重新验证可用绑定和重载。[本批验证](#p1-mpi-绑定计数及分块) | 按部署环境记录 1/2/4 ranks 正常与失败路径、空分区、分块结果和计数范围。 |
| `grasp_shared` · P1 | LP64 配置、稀疏矩阵存储、活动 SCF 系数和基组数组所有权已进行本地验收。`alloc` 约定传入无活跃目标的指针，先 `nullify` 再 `associated` 无法检测误传的活目标。按调用方选择拥有数据的 `allocatable` 或有明确生存期的借用视图。[LP64 配置验收](#p1-blaslapack-lp64-接口与配置验收) · [稀疏矩阵验收](#p1--稀疏矩阵存储与数组所有权) · [分配接口](../../src/lib/grasp_shared/array_memory_management.f90#L11) | 对改动的分配路径测试零长度、非 1 下界、重分配保留重叠区、释放后再用。维度扩大须在算术前转换，并在转入默认整数前检查上界。 |
| `core` · P1，写入器多实例 P2 | 混合系数写入器约定每进程一次写一个流，MPI 由 root 写。已核对打开、写块、刷新、完成标记及关闭路径，并覆盖短写和根进程失败。双正交代码在文件末尾关闭输入和输出。有交错输出的真实需求时再引入显式写入器对象。径向网格与稀疏矩阵乘则按已确认的共享写入依赖收紧接口。[混合流](../../src/lib/core/mixing_file_format.f90#L25) · [双正交消费者](../../src/apps/rbiotransform/read_biorthogonal_mcp.f90#L188) · [网格](../../src/lib/core/radial_grid_setup.f90) | 正常文件可由独立解析器读回。核对块目录、基组摘要、根顺序和完成标记，验证中断或截断文件及 root 写入失败。保留各程序活动网格、填充点与精度初始化时序，并用同输入核对径向数组与矩阵乘结果。 |
| `csf_descriptor` 与 Fortran reader · P1 | 现行 V3 跨语言接口已有版本、容量及分块检查。`descriptor_source` 拥有 Rust 句柄，调用方须显式关闭并避免复制活句柄。已覆盖不同基组交替读取、句柄关闭及失败后的重开和替换。[Fortran 所有权](../../src/lib/core/csf_descriptor_reader.f90#L19) · [C 绑定](../../src/lib/core/csf_descriptor_c_bindings.f90) · [Rust FFI](../../src/lib/csf_descriptor/src/descriptor_ffi.rs) | 双基组交替读取、失败后重开、重复关闭及错误摘要均可诊断。每个所有者只释放一次。块顺序、占据与耦合码逐项一致，分块读取记录峰值内存。 |
| `davidson_eigensolver` · P1 | 保留显式矩阵乘回调，核对内部 assumed-size 数组的实际形状、下界和调用方式，先处理影响索引校验的入口。已用近简并双根夹具验证重启、返回向量重用和按根的残差判定。[回调](../../src/lib/davidson_eigensolver/davidson_eigensolver.f90) · [稀疏初值](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90) | 扩展已有解析矩阵和提前收敛用例，增加近简并多根、重启向量、空分区。分别比较每根本征残差、正交性和近简并子空间。 |
| `spin_angular` · P1 核查，P2 按需整理 | 跟踪一个 bra/ket 配对的可变工作区。Wigner 6j 的显式缓存已有线程私有语义，检查冷启动、复用和线程数变化。只读表与可变缓存分别处理。[配对数据](../../src/lib/grasp_shared/csf_pair_angular_state.f90) · [6j 缓存](../../src/lib/spin_angular/wigner_6j.f90#L33) | 保留 NMTEJJ 整数对照，覆盖高 J、交换相位和解析选择定则零值。将角系数逐对读回，在实际启用 OpenMP 的消费者上另做多线程对照。 |
| `breit_angular_support` · P1 | 明确 Breit 系数累加与回调的输入输出。非法轨道组合和零张量秩已有可诊断的串行失败和 MPI 终止回归。[径向系数](../../src/lib/breit_angular_support/breit_radial_integral_coefficients.f90#L311) · [回调回归](../../test/f2018/breit_procedure_hooks.f90) | 对照直接项、交换项、开关前后矩阵元素和相位；非法组合产生可识别的失败。用固定轨道 RCI 确认库测试覆盖的分支可在生产路径触发。 |

跨库审查应随涉及的调用链执行。以实际编译目标、宏和模块目录重新收集编译器及 Fortitude/fortls 诊断，优先逐条查明未初始化、越界、窄化、接口与跨调用局部状态。改动 `INTENT`、可见性、固定单元号或 `USE ONLY` 前检查关键字调用、重导出、回调和错误清理。早期提示数和 529 项命名映射属于各自快照。GNU 16.2.1/Linux 历史复采的 37 个新增诊断身份中，9 个是 `-Wmaybe-uninitialized`，复核时记录赋值路径和可达输入。[告警基线](../reference/f2018/evidence-0019-gnu16-warning-baseline.md) · [静态审查及工具限制](../reviews/evidence-0012-fortls-fortitude-audit.md) · [命名映射](../reference/non_csfg_variable_renames_2026-11-16.json)

## 应用程序实施与退出条件

同一行的串行、MPI 或 phase 变体分别登记对照结果。各变体先与自身版本基线比较，再解释跨变体的归约顺序、存储方式和计算分支差异。

| 应用目标 | 工作项 | 退出条件 |
| --- | --- | --- |
| `rwfnestimate` · P1 | 固定点核、有限核和自定义网格，分别核对导入轨道、Thomas–Fermi、氢样初猜及正交化路径。复用既有 Ni/Zn 高角动量夹具。[入口](../../src/apps/rwfnestimate/rwfnestimate_main.f90) | 比较轨道身份、节点数、同 κ 正交、归一化及相位对齐后的 P/Q，输出可供首轮 RMCDHF 读入。坏轨道元数据、点数与表长有明确诊断。 |
| `rangular`、`rangular_mpi` · P1 | 核对配对工作区、计数、排序、续写和空分区。串行按 MCP 通道读取，MPI 按 rank 使用索引分片。[生成路径](../implemented/design-0017-rangular-flow-and-state.md) · [MPI MCP](../../src/lib/mpi/mpi_indexed_mcp.f90) | 逐对称块、CSF 对和积分标签核对稀疏位置、系数与相位。覆盖首阶限制、串行复用开关、错基组和分片不匹配，由相应 RMCDHF 消费端读回。 |
| `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` · P0 首轮夹具，P1 多根 | 分别冻结 EOL 根、权重、初始轨道、变化轨道、阻尼和停止参数。跟踪每块矩阵装配、每轮混合更新、轨道更新次序以及缓存或收缩失效条件。普通 MPI 按当前核势入口验收，另三种目标覆盖支持的 NP 关闭、零势和非零势。[SCF 状态及顺序](../implemented/design-0019-rmcdhf-flow-and-state.md) · [变体差异](../reviews/evidence-0007-rmcdhf-memory-variant-comparison.md) · [核极化范围](../implemented/design-0022-rmcdhf-nuclear-polarization.md) | 每个变体保存同输入下的固定轮次轨迹与独立收敛轨迹。逐轮核对能量、根与权重、系数、残差、P/Q 和停止原因。MPI 分别生成与本次 rank 布局匹配的 MCP；记录各 rank 峰值内存，检查 memory 路径是否因所有权调整引入多余副本。 |
| `rci`、`rci_mpi` · P0 续算，P1 固定轨道 | 首批完成读取及局部失败回归。随后以同轨道、同基组、同算符，分别覆盖小稠密块、稀疏矩阵乘与磁盘续算，逐一切换 Coulomb、Breit 和实际支持的质量位移、QED 选项。[块求解](../../src/apps/rci/solve_hamiltonian_blocks.f90) · [已有 Breit 回归](../../test/rci_breit.py) | 核对块维度、矩阵元素和对称性、每根残差、总能量及激发能；冷启动与续算对照态身份与系数。保存最大差异的块、根和元素位置。 |
| `rbiotransform`、`rbiotransform_mpi` · P1 | 主程序调用的 `prepare_biorthogonal_subshells` 已判定正序、逆序、两侧不一致及非法排列。覆盖这一实现的退出语义，并核查双基组映射、重叠矩阵、临时工作区和角缓存身份。[实际顺序判断](../../src/apps/rbiotransform/prepare_biorthogonal_subshells.f90#L180) · [MPI 入口](../../src/apps/rbiotransform_mpi/rbiotransform_mpi_main.f90#L170) | 正序与共同逆序各有正例，两侧不一致和混乱排列各有负例。用已知双正交关系验证交叉重叠和相应矩阵元，覆盖轨道数不同、重叠矩阵奇异及错误缓存。核对输出流的来源基组摘要。 |
| `rtransition`、`rtransition_mpi`、`rtransition_phase` · P1 验收，P2 按需去重 | 固定初末态映射、相位约定、跃迁能量和角缓存。对 phase 版单独跟踪带符号矩阵元写出。公共计算去重以各分支已有对照为前提。[普通入口](../../src/apps/rtransition/README.md) · [带符号输出](../../src/apps/rtransition_phase/README.md) | 完成双正交到跃迁的允许、禁戒和弱线用例。分别比较电多极两规范、磁多极的矩阵元和线强；固定频率比较率，再检验使用计算能量的频率因子。phase 版对齐初末态相位后核对符号。 |
| `gj`、`gj_mpi` · P1 | 对 MCDHF 与 CI 输入匹配相同 ASF，固定算符与修正选项，沿各入口核对矩阵元和归约。[g 因子计算](../../src/apps/gj/compute_gj_only.f90) | 小模型对照解析 LS 极限，真实夹具逐态对照 Landé 因子、J 与宇称。记录串行及 MPI 的差异。 |
| `rhfs`、`rhfs_mpi` · P1 | 跟踪约化矩阵元工作区的分配与释放，固定同位素、核自旋、核矩和单位。[矩阵元实现](../../src/apps/rhfs/calculate_hyperfine_matrix_elements.f90) | 各态分别对照 A、B 和相位一致的非对角元。覆盖 CI 与非 CI 路径，检查 MPI 归约结果及输出单位。 |
| `ris4`、`ris4_mpi`、`fical` · P1 | RIS 缓存已有缺失、截断和重用回归，补同维度错误基组或选项的核查。FICAL 对文本输入分段检查读状态与量纲。[缓存读入](../../src/apps/ris4/open_isotope_shift_angular_cache.f90) · [两态位移入口](../../src/apps/ris4/isotope_shift_transition_main.f90) | 分别核对 NMS、SMS、场位移及其总和；缓存重算和复用结果一致。FICAL 比较两态差的符号、单位和截断记录的报错。 |
| `rdensity`、`rdensity_mpi` · P1 | 检查密度矩阵形状、MPI 归约的 `size(matrix)` 计数和自然轨道求解时的工作区。[MPI 归约](../../src/apps/rdensity/density_parallel_mpi.f90#L55) | 验证密度矩阵对称性、迹与电子数、占据数、自然轨道正交及流文件读回。比较空 rank 和简并子空间，沿用现有串行及 1/2/4-rank 密度用例。 |
| `jj2lsj_omp` · P1 验收 | 分开核对控制器状态、只读表和线程私有缓存。关闭 OpenMP 时同一目标使用串行后端。[目标定义](../../src/apps/jj2lsj/CMakeLists.txt) · [已有 Gd I 回归](../reviews/evidence-0002-jj2lsj-gd-i-validation.md) | 串行与 OpenMP 1/2/4 线程对照完整系数、排序、标签和提前停止路径；在改动缓存生存期时加重复调用或线程组变化用例。 |
| `rlevels` · P2 | 复用混合流 reader，核对多文件能级排序、可选 LSJ 标签匹配和高 J 输出。[能级入口](../../src/apps/rlevels/rlevels.f90) | 逐态核查绝对能量、激发能、相邻能级间隔与标签归属。覆盖同能不同对称性、缺标签和损坏记录，结构化输出继续通过现有文件回归。 |

## 数值对照协议

1. **固定对照条件。** 选取可重放的小型单根模型、含近简并块的多根模型，以及一套来源可追溯的真实多根输入，例如 Zn I 奇宇称小 MR。记录核模型、网格与精度、CSF 及块顺序、初猜、变分轨道、每块根选择、EOL 权重和算符开关。保存提交及构建选项、编译器/MPI/BLAS 版本、rank 和线程布局、输入摘要。近简并簇和容差在参考构建上预先确定。
2. **先比身份，再比数值。** 块的 J 与宇称、CSF 顺序、轨道量子数和离散映射精确比对。非简并态在相同对称块内按重叠匹配，比较向量和带符号矩阵元前对齐相位。近简并簇比较投影子空间，孤立的根编号或能量排序不足以确认态身份。轨道比较包含 P/Q 两分量、节点、归一化和同 κ 正交。
3. **为不同问题保留各自的判据。** 下表摘取历史方案中条件良好小夹具的候选阈值。实施时先用足够精度的输出检查基准稳定性，记录最终采用的阈值和依据。生产模型精度、活性空间收敛、I/O 性能另行评估。[候选阈值与单位](../archive/plan-0002-f2018-migration-gnu15.md#L223)

| 量 | 判读方法 |
| --- | --- |
| RCI 总能量、激发能 | 固定轨道小夹具分别比较两者的绝对差，各自候选上限为 `1e-10 Ha`。精细结构间隔另列。 |
| CI 本征方程残差 | 对每个根计算 `norm(Hc − Ec)₂ / ((norm(H)F + abs(E)) × norm(c)₂)`，候选上限为 `1e-10`。分母为零单独处理。 |
| CI 系数归一化与正交 | 计算 `max(abs(CᵀC − I))`，候选上限为 `1e-12`。同时保留最差的根或根对。 |
| SCF 固定轮次和正常收敛 | 固定轮次比较相同迭代位置的能量、系数、轨道和分支。收敛用例要求两边各自满足冻结的停止条件，记录停因、轮次和最终态匹配。达到迭代上限仅用于固定轮次结论。 |
| 性质和跃迁 | 固定核矩、算符和单位，分项对照。弱线、解析零或强抵消项使用绝对误差；电跃迁的长度及速度形式分别与自身基线比较。 |

## 执行入口与验收记录

以下示例从 Suite 仓库根目录运行。先在目标机加载匹配的 GNU Fortran 16.2、Open MPI 和 BLAS/LAPACK；示例使用工作区共享 Python 环境的 V3 `rcsfs`。每次建立独立构建目录，临时文件集中在工作区的 `temp/tests`，以满足双基组 MPI 流程的路径长度限制。第一组命令运行**已有**接口及工作流回归。新增的 RCI 坏记录、单分片失败和首轮实际二进制对照应在实现后加入对应测试组。

```fish
set -gx FC (command -v gfortran)
set -gx OMPI_FC "$FC"
$FC -dumpfullversion
mpifort -dumpfullversion
set suite_venv (realpath ../graspkit-tools/.venv)
set suite_python "$suite_venv/bin/python"
set -gx RCSFS_TEST_CLI "$suite_venv/bin/rcsfs"
mkdir -p ../temp/tests
set suite_run (realpath (mktemp -d ../temp/tests/suite-f2018.XXXXXXXX))
set -gx TMPDIR (realpath ../temp/tests)
./configure.sh --debug --profile f2018-check --build-dir "$suite_run/check" -DRCSFS_TEST_CLI="$RCSFS_TEST_CLI"
cmake --build "$suite_run/check" --parallel 4
ctest --test-dir "$suite_run/check" --no-tests=error --output-on-failure -R '^(f2018_linear_algebra_contracts|f2018_array_memory|f2018_davidson_contracts|state_workflow_(serial|mpi))$'
```

跨共享接口、绑定或发布候选的完整配置验收使用同一环境和以下命令。矩阵覆盖三种 profile，各有 serial、MPI 和 jj2lsj OpenMP 的八种组合。保留 `matrix-validate.json`、每个配置的 CTest 结果和严格配置的逐文件告警比对。某个组合单独通过时，按该组合报告。[矩阵执行器](../../test/f2018/matrix.py) · [警告预算与测试范围](../../test/f2018/README.md)

```fish
"$suite_python" test/f2018/matrix.py "$suite_run/matrix" \
    --rcsfs-cli "$RCSFS_TEST_CLI" --parallel-configs 2 --jobs 3
```

每批在结果记录中写明基准与候选提交、受影响的程序、实际运行的配置、输入与日志位置、通过/失败/跳过项及未覆盖路径。数值对照同时报告最差差异的块、态、单位和终止原因。对已确认的错误保存修复前的表现、修复后的正反例和合法输入对照；对仅由诊断提示的风险保存调用链及可达性判断。完成相应验收后再更新本计划的状态。

| 交付批次 | 完成时应能查到的结果 | 编制时状态 |
| --- | --- | --- |
| 首批 · 基线与续算 | V3 实际首轮对照、串行/MPI 合法续算、坏记录和单 rank 失败报告 | GNU Fortran 16.2.1/OpenBLAS/Check 下，两套独立构建的 `rmcdhf`、`rmcdhf_mem` 首轮科学输出一致。串行和 Open MPI 5.0.11 的 1/2/4-rank 续算回归已完成。 |
| MPI、ABI 与共享库 | 绑定和计数用例、链接库身份、错误传播、所有权及共享核差异报告 | 本机已增加 32 位边界、CSF 基组及活动系数所有权、混合流故障、Davidson 近简并和 Breit 回归。其余共享接口与部署环境继续按批核验。 |
| 各程序的科学数值回归 | 四种 RMCDHF 各自的固定轮次/收敛记录、固定轨道 RCI、双正交及性质与跃迁分项对照 | 四个 RMCDHF 变体已有小型多根对照，固定轨道 RCI 已核对小稠密块残差。生产规模多根、RCI 存储分支及性质与跃迁完整对照仍待验收。 |
| 整体验收 | 候选提交的 24 配置结果、正确性告警逐项处理结论及目标部署环境报告 | 历史矩阵供追溯，候选版本另行签收 |
