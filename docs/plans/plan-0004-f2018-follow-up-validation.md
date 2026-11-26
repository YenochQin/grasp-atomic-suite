# Suite F2018 现代化实施与验收计划

**状态：待实施。** 2026-11-26 按 Suite `3.4dev1` 的 `f211a0a35a957fab85c65a8a514d779183d0d59e` 复核。第一批处理 RCI 串行与 MPI 的续算读取、单 rank 失败传播，并建立能从 V3 生产入口运行的版本对照。后续批次处理 MPI 绑定、共享数值接口和多根计算链。本次编制只核对了源码、测试定义与历史记录，表中验收项均需由相应实施批次给出结果。

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

**现状。** 独立首轮比较器在新目录只写 `isodata` 和 `rcsf.inp`，随后启动 `rwfnestimate → rangular → RMCDHF`。它还严格核对输出文件清单。当前初猜从 header 取轨道表，角系数和 SCF 从 V3 文件对取基组。CTest 中的 `f2018_first_sweep_driver` 测试比较工具本身，单独跑它不能取得两个实际构建的首轮结果。[比较器](../../test/f2018/compare_first_sweep.py#L148) · [CTest 注册](../../test/f2018/CMakeLists.txt#L83) · [V3 夹具生成器](../../test/native_fixture.py)

**实施。** 用现有夹具生成器在独立目录生成配对的描述符和 header，给两侧提供相同的输入字节、状态名和交互应答。更新文件清单，既识别生成的输入，也检查输出完整。用相同编译器、数学库、优化选项及线程数分别构建基准和候选，存下输入、可执行文件摘要、标准输出和错误输出。先用当前源码的两个独立构建验证比较器能稳定识别一致结果，再以坏签名、缺块或改动一个科学值验证它能报告差异。历史版需要旧格式时，单独固定参考输入转换和字段映射。

**交付与退出。** 对 `rmcdhf` 和 `rmcdhf_mem` 分别取得可重放的一轮结果。核对 MCP 块与记录、轨道身份及 P/Q、ASF 根、能量和系数。明确预期停因为迭代上限。较大的真实多根算例在涉及相应内核改动前冻结，首轮软件回归可以独立完成。[现有工作流夹具](../../test/test_state_workflow.py) · [首轮诊断](../../test/f2018/compare_first_sweep.py#L110)

### P0 · 修正 RCI 续算记录的读取边界

**源码证据。** `rci` 与 `rci_mpi` 的 `read_ci_restart` 在一条顺序无格式 `READ` 中先取得 `radial_point_count`，随即以它为上界读入核势数组，之后才检查 `max_radial_points`。核模型参数个数也直接决定同一记录中数组的读入长度，而数组容量为 2。每轨道尾点从前一记录读出后，下一次读取直接用作 P/Q 上界。源码已能确定校验顺序，越界或异常退出在具体构建上的表现待复现。[串行读取](../../src/apps/rci/read_ci_restart.f90#L47) · [MPI 读取](../../src/apps/rci_mpi/read_ci_restart.f90#L50) · [核参数容量](../../src/lib/grasp_shared/nuclear_model_parameters_state.f90) · [径向容量](../../src/lib/grasp_shared/cal_constants.f90)

1. 从写端列出每条记录的字段、计数、目标数组和后继使用者。用实际写端生成的小型合法文件作为正例。点核的参数个数为 0，现有有限核模型为 2。核对活动网格、含 10 个填充点的工作容量、每轨道尾点，以及 CSF 和对称块计数。[模型输入](../../src/lib/core/isotope_data_records.f90) · [串行写端](../../src/apps/rci/read_ci_options.f90#L215)
2. 在依赖计数的读取前验证取值和目标形状。对“计数与数组同记录”的情况，可先仅读计数，检查成功后 `BACKSPACE` 回到该记录再完整读取，对探读、定位和重读分别检查 `iostat`。顺序无格式文件中两次普通 `READ` 会消费两条记录，须用既有合法文件验明记录位置。尽量将记录解析供串行和 MPI 共用，由入口决定退出方式。
3. 按文件名、记录及适用的块或轨道报告首个错误。用合法边界值、越界计数、负计数、截断记录和非法轨道尾点做回归。正例做同模型冷启动与续算的态、能量和系数对照。负例要求程序主动给出可定位的非零失败，避免仅依赖 Debug 的运行时检查。
4. 单列基组身份问题。现有头记录可核对电子、CSF、轨道和块数，同维度文件仍可能对应不同的 CSF 排列。先追踪完整写入、读取与分片信息，再评估是否增加摘要或版本。若需改动磁盘协议，另列读写端覆盖、既有文件处置及同维异基组回归，把该项与读取边界修复独立提交。

### P0 · 处理 RCI 的单 rank 分片失败

**源码证据。** `rci_mpi` 通过 root 广播续算选择，各 rank 随后打开自己的 `rciNNN.res`。某一片缺失、打开失败或签名错误时，该 rank 可单独进入 `finalize_mpi_runtime`。这个过程先对 `MPI_COMM_WORLD` 做 barrier，再 finalize。其他 rank 的控制流可以继续读各自分片，调用序列存在分歧。实际是否挂起、报出什么以及启动器退出码，须在指定 MPI 上测量。[分片入口](../../src/apps/rci_mpi/connect_ci_restart_file.f90#L50) · [退出实现](../../src/lib/mpi/mpi_runtime.f90#L54)

**实施。** 在所有 rank 均能到达的阶段汇总本地状态，或在不可恢复的本地 I/O 错误上复用已有 `mpi_io_errors` 终止通路。保证出错 rank 能报告文件、rank 和阶段。沿 RCI 的头记录与矩阵记录调用链排查其他依赖集体到达的错误出口。正常完成仍经应用的集体收尾路径。[MPI I/O 错误处理](../../src/lib/mpi/mpi_io_errors.f90) · [续算矩阵乘](../../src/apps/rci_mpi/multiply_restart_hamiltonian.f90)

**交付与退出。** 以 2 ranks 先生成匹配的文件，分别只让 root 或一个 worker 的分片缺失、签名错误或记录截断。调用真正的续算入口，检查启动器非零退出、输出命中预期分片和阶段、所有进程在测试超时前结束。将“命令超时”和“收到预期失败”记录为两种不同结果。用 1/2/4 ranks 的合法分片覆盖正常退出及结果读回，每次按对应 rank 布局准备文件。P0 的实现可以使用当前 MPI 绑定。

## 后续批次与共享库

先用第一批夹具保护读写和失败语义，再迁移通信接口。改动共用矩阵、轨道或角系数之前，先取得相应消费者的固定输入数值基线。单个子项目可按自身依赖推进。以下路径以实现库为单位；描述符的 Rust 端和 Fortran 调用端共同验收。

| 位置与优先级 | 实施范围和已有基础 | 验收证据 |
| --- | --- | --- |
| `mpi` · P1 | `mpi_state` 目前转发 `use mpi` 并保存整数状态数组。先验证部署 MPI 的 `mpi_f08` 和所需重载，再按通信调用链迁移 communicator、datatype、status 及实际使用的其他句柄。常规 rank 仍是整数。分别处理元素计数、地址位移和文件偏移。CSF 广播已有有界分块，可复用于有相同约束的路径。[状态模块](../../src/lib/mpi/mpi_state.f90) · [分块实现](../../src/lib/mpi/mpi_csf_collectives.f90#L68) · [MPI 绑定说明](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node465.htm) | 正常与错误路径均运行 1/2/4 ranks 中适用的组合。覆盖空分区、限小分块、非连续数据的打包或已验证的传递方式、整数边界和窄化拒绝。在支持大计数的部署上验证对应重载；在有界实现上测试跨分块结果与原顺序。 |
| `grasp_shared` · P1 | 现有 BLAS/LAPACK 声明及实际运算测试已按 LP64 建立。补齐 CMake 对整数宽度的选择和部署验证。沿稀疏矩阵列尾偏移、SCF 系数与基组指针追踪分配者和有效期。`alloc` 约定传入无活跃目标的指针，先 `nullify` 再 `associated` 无法检测误传的活目标。先核对调用方，再选择拥有数据的 `allocatable` 或有明确生存期的借用视图。[BLAS 回归](../../test/f2018/linear_algebra_contracts.f90) · [分配接口](../../src/lib/grasp_shared/array_memory_management.f90#L11) · [稀疏存储](../../src/lib/grasp_shared/sparse_hamiltonian_storage.f90) | 验证实际链接库与整数、浮点 ABI 相符，记录不匹配环境的拒绝方式。对改动的分配路径测试零长度、非 1 下界、重分配保留重叠区、释放后再用。维度扩大须在算术前转换，并在转入默认整数前检查上界。 |
| `core` · P1，写入器多实例 P2 | 混合系数写入器约定每进程一次写一个流，MPI 由 root 写。当前优先追踪打开、写块、刷新、完成标记及关闭失败能否传回正确入口。双正交代码在文件末尾关闭输入和输出。有交错输出的真实需求时再引入显式写入器对象。径向网格与稀疏矩阵乘则按已确认的共享写入依赖收紧接口。[混合流](../../src/lib/core/mixing_file_format.f90#L25) · [双正交消费者](../../src/apps/rbiotransform/read_biorthogonal_mcp.f90#L188) · [网格](../../src/lib/core/radial_grid_setup.f90) | 正常文件可由独立解析器读回。核对块目录、基组摘要、根顺序和完成标记，验证中断或截断文件及 root 写入失败。保留各程序活动网格、填充点与精度初始化时序，并用同输入核对径向数组与矩阵乘结果。 |
| `csf_descriptor` 与 Fortran reader · P1 | 现行 V3 跨语言接口已有版本、容量及分块检查。`descriptor_source` 拥有 Rust 句柄，调用方须显式关闭并避免复制活句柄。核查双基组装载时句柄、摘要和当前块是否始终对应正确的一侧。[Fortran 所有权](../../src/lib/core/csf_descriptor_reader.f90#L19) · [C 绑定](../../src/lib/core/csf_descriptor_c_bindings.f90) · [Rust FFI](../../src/lib/csf_descriptor/src/descriptor_ffi.rs) | 双基组交替读取、失败后重开、重复关闭及错误摘要均可诊断。每个所有者只释放一次。块顺序、占据与耦合码逐项一致，分块读取记录峰值内存。 |
| `davidson_eigensolver` · P1 | 保留显式矩阵乘回调，核对内部 assumed-size 数组的实际形状、下界和调用方式，先处理影响索引校验的入口。检查工作区重用和每根的收敛标记。[回调](../../src/lib/davidson_eigensolver/davidson_eigensolver.f90) · [稀疏初值](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90) | 扩展已有解析矩阵和提前收敛用例，增加近简并多根、重启向量、空分区。分别比较每根本征残差、正交性和近简并子空间。 |
| `spin_angular` · P1 核查，P2 按需整理 | 跟踪一个 bra/ket 配对的可变工作区。Wigner 6j 的显式缓存已有线程私有语义，检查冷启动、复用和线程数变化。只读表与可变缓存分别处理。[配对数据](../../src/lib/grasp_shared/csf_pair_angular_state.f90) · [6j 缓存](../../src/lib/spin_angular/wigner_6j.f90#L33) | 保留 NMTEJJ 整数对照，覆盖高 J、交换相位和解析选择定则零值。将角系数逐对读回，在实际启用 OpenMP 的消费者上另做多线程对照。 |
| `breit_angular_support` · P1 | 明确 Breit 系数累加与回调的输入输出。跟踪非法轨道组合触发的库内 `STOP`，由相应串行或 MPI 入口负责故障终止。[径向系数](../../src/lib/breit_angular_support/breit_radial_integral_coefficients.f90#L311) · [回调回归](../../test/f2018/breit_procedure_hooks.f90) | 对照直接项、交换项、开关前后矩阵元素和相位；非法组合产生可识别的失败。用固定轨道 RCI 确认库测试覆盖的分支可在生产路径触发。 |

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

以下示例从 Suite 仓库根目录运行。先在目标机加载匹配的 GNU Fortran 16.2、Open MPI 和 BLAS/LAPACK；示例使用工作区共享 Python 环境的 V3 `rcsfs`。每次建立独立目录，供后续批次保留构建树、临时文件和结果。第一组命令运行**已有**接口及工作流回归。新增的 RCI 坏记录、单分片失败和首轮实际二进制对照应在实现后加入对应测试组。

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
mkdir -p "$suite_run/tmp"
set -gx TMPDIR "$suite_run/tmp"
./configure.sh --debug --profile f2018-check --build-dir "$suite_run/check" -DRCSFS_TEST_CLI="$RCSFS_TEST_CLI"
cmake --build "$suite_run/check" --parallel 4
ctest --test-dir "$suite_run/check" --no-tests=error --output-on-failure -R '^(f2018_linear_algebra_contracts|f2018_davidson_contracts|state_workflow_(serial|mpi))$'
```

跨共享接口、绑定或发布候选的完整配置验收使用同一环境和以下命令。矩阵覆盖三种 profile，各有 serial、MPI 和 jj2lsj OpenMP 的八种组合。保留 `matrix-validate.json`、每个配置的 CTest 结果和严格配置的逐文件告警比对。某个组合单独通过时，按该组合报告。[矩阵执行器](../../test/f2018/matrix.py) · [警告预算与测试范围](../../test/f2018/README.md)

```fish
"$suite_python" test/f2018/matrix.py "$suite_run/matrix" \
    --rcsfs-cli "$RCSFS_TEST_CLI" --parallel-configs 2 --jobs 3
```

每批在结果记录中写明基准与候选提交、受影响的程序、实际运行的配置、输入与日志位置、通过/失败/跳过项及未覆盖路径。数值对照同时报告最差差异的块、态、单位和终止原因。对已确认的错误保存修复前的表现、修复后的正反例和合法输入对照；对仅由诊断提示的风险保存调用链及可达性判断。完成相应验收后再更新本计划的状态。

| 交付批次 | 完成时应能查到的结果 | 编制时状态 |
| --- | --- | --- |
| 首批 · 基线与续算 | V3 实际首轮对照、串行/MPI 合法续算、坏记录和单 rank 失败报告 | 已完成源码复核，修复与运行待实施 |
| MPI、ABI 与共享库 | 绑定和计数用例、链接库身份、错误传播、所有权及共享核差异报告 | 部分基础回归已注册，迁移及部署验收待实施 |
| 各程序的科学数值回归 | 四种 RMCDHF 各自的固定轮次/收敛记录、固定轨道 RCI、双正交及性质与跃迁分项对照 | 有既有夹具和专项记录，完整多根对照待实施 |
| 整体验收 | 候选提交的 24 配置结果、正确性告警逐项处理结论及目标部署环境报告 | 历史矩阵供追溯，候选版本另行签收 |
