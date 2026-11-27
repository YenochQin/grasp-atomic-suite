# Suite Fortran 命名、注释与 README 审查

2026-11-22，源码基准为 `9d0a85fd`，审查和修复均基于当时的工作树。1,002 个生产 Fortran 文件均符合小写 `snake_case` 的文件名格式。下列各项已修正，修改限于文件名、应用内部标识符、注释、文档及对应的源码快照记录。

## 检查范围

修复前的文件名扫描覆盖 `src/apps/` 的 698 个文件和 `src/lib/` 的 304 个文件。`grasp_shared` 的 110 个 `.f90` 文件均定义同名模块。20 份应用 README 覆盖全部 698 个应用源码文件。35 份 README 中扫描到的 1,262 个仓库内链接均可定位。标识符的含义、注释和计算流程说明按实现抽查。

## 文件名和标识符

- `rbiotransform` 和 `rbiotransform_mpi` 的原 `check_orbital_order.f90` 文件名暗示会检查轨道顺序。两份文件均改为 `check_orbital_order_stub.f90`，保留历史模块名、过程接口和直接返回的行为。应用 README 与迁移记录已更新路径，实际轨道顺序检查仍由 `prepare_biorthogonal_subshells` 执行。[串行占位过程](../../src/apps/biorthtransform/check_orbital_order_stub.f90) · [MPI 占位过程](../../src/apps/biorthtransform_mpi/check_orbital_order_stub.f90) · [轨道顺序处理](../../src/apps/biorthtransform/prepare_biorthogonal_subshells.f90#L189)
- 串行和 MPI 的 `tabulate_slater_integrals` 将原 `n` 和 `key` 改为 `slater_integral_count` 和 `orbital_index_base`。基数仍为 `orbital_count + 1`。串行主程序和 MPI 归约包装过程中的计数名称也已同步。[串行 Slater 积分表](../../src/apps/rci/tabulate_slater_integrals.f90#L27) · [MPI Slater 积分表](../../src/apps/rci_mpi/tabulate_slater_integrals.f90#L27) · [MPI 归约](../../src/apps/rci_mpi/tabulate_slater_integrals_mpi.f90#L27)
- `rwfnestimate` 将原 `b1` 至 `b4` 改为表示当前点及沿积分方向前三个已知点的权重名称。声明旁说明这些权重作用于径向大、小分量采样值，导数权重来自 `starting_step_weights`。[径向传播常量](../../src/apps/rwfnestimate/integrate_homogeneous_radial_segment.f90#L34)

## 源码注释

- 两份轨道顺序占位过程的文件头现在明确说明入口直接返回，并移除了宣称会执行检查的旧注释。[串行文件头](../../src/apps/biorthtransform/check_orbital_order_stub.f90#L1) · [MPI 文件头](../../src/apps/biorthtransform_mpi/check_orbital_order_stub.f90#L1)
- 四种普通 RMCDHF 的 SCF 循环将旧的 `.rwf` 注释改为说明每轮写入调用方选定的输出路径。命名计算使用 `name.wfn`，空状态名使用 `rwfn.out`。四份导入源码的注释快照已同步。[串行循环](../../src/apps/rmcdhf/run_scf_iterations.f90#L203) · [串行内存版循环](../../src/apps/rmcdhf_mem/run_scf_iterations.f90#L203) · [MPI 循环](../../src/apps/rmcdhf_mpi/run_scf_iterations_mpi.f90#L226) · [MPI 内存版循环](../../src/apps/rmcdhf_mem_mpi/run_scf_iterations_mpi.f90#L226) · [状态路径](../../src/lib/core/calculation_state_files.f90#L32)

## README 与实现的差异

### FICAL 编译示例

旧的单文件编译命令在干净目录报找不到 `calculation_state_files.mod`。RIS4 README 已改用配置后的 CMake `fical` 目标。该目标链接 `core`，提供 FICAL 读取上态和下态名称所需的模块及实现。[RIS4 README](../../src/apps/isotope/README.md#L26) · [FICAL 入口](../../src/apps/isotope/isotope_shift_transition_main.f90#L4) · [CMake 目标](../../src/apps/isotope/CMakeLists.txt#L1)

### gj 示例输入

按旧示例直接运行 `gj` 会因缺少 `test_header.toml` 退出。中英文 README 已加入可复现步骤，用支持 V3 的 `rcsfs` 将 `test.c` 转为 `test_desc.parquet` 和配套头文件，然后用归档的 `test.m` 和 `test.w` 运行 MCDHF 模式。通用 `gj_mpi` 命令同时说明工作目录、可执行文件路径及 `--ci` 需要的 `name.cimix` 或 `name.cm`。[英文示例和样本](../../README.md#L217) · [中文示例和样本](../../README_ZH.md#L194) · [状态文件解析](../../src/lib/core/calculation_state_files.f90#L113) · [参考输出](../../data/test.gj)

### CSF 编码和共享库过程名

- Core README 已改用 V3 字段顺序和编码。总 `2J+1` 与宇称分别写入 `csf_total_j_slot` 和 `csf_parity_slot`，查询函数直接读取各槽位。原先的带符号合并编码描述已替换。[Core README](../../src/lib/core/readme.md#L349) · [V3 字段定义](../../src/lib/csf_descriptor/src/descriptor_transport.rs#L19) · [Fortran 安装槽位](../../src/lib/core/csf_descriptor_reader.f90#L280) · [查询函数](../../src/lib/core/csf_quantum_queries.f90#L76)
- Core、Davidson 和自旋角库 README 已更新旧过程名，包括 `radgrd` 对应的 `initialize_radial_grid`、`mgs_nrm` 对应的 `orthonormalize_davidson_vectors`，以及 `onescalar`、`oneparticlejj` 对应的两个单电子系数入口。共享库 README 中的旧 Slater 制表名称也已更新。[Core README](../../src/lib/core/readme.md#L9) · [Davidson README](../../src/lib/davidson_eigensolver/readme.md#L5) · [自旋角库 README](../../src/lib/spin_angular/readme.md#L5) · [共享库 README](../../src/lib/grasp_shared/readme.md#L9)

### MPI 内存版 RMCDHF 的停止条件和输出

应用 README 已说明 rank 0 每轮写径向波函数，轨道自洽和相对加权能量变化可分别触发停止。达到最大迭代次数时也会结束并打印提示，收敛状态应结合迭代报告判断。[MPI 内存版 README](../../src/apps/rmcdhf_mem_mpi/README.md#L13) · [SCF 迭代与停止判断](../../src/apps/rmcdhf_mem_mpi/run_scf_iterations_mpi.f90#L190)

### 顶层贡献指南的文件命名

中英文 README 已将旧的 `*_I.f90` 和 `*_C.f90` 要求更新为现行约定。应用按职责使用小写 `snake_case` 文件名，迁移后的模块可保留 `*_interface` 名称，`grasp_shared` 模块与文件同名。[英文指南](../../README.md#contributor-notes) · [中文指南](../../README_ZH.md#开发注意事项) · [应用命名规范](../../AGENTS.md#coding-style--naming-conventions) · [共享模块约定](../../AGENTS.md#shared-libraries-and-retained-programs)

## 验证记录

修复前，在 macOS 上用 GNU Fortran 16.2.0 构建 `gj`、`gj_mpi`、`rlevels` 和 `fical`，选跑的 6 项 CTest 全部通过。旧档直接运行 `gj` 和 FICAL 单文件编译的失败均已复现。

修复后，在同一编译器上完成 CMake Release 全量构建，包含四种普通 RMCDHF、串行和 MPI RCI、`rwfnestimate`、双正交程序及相关测试目标。

- 对 11 份修改过的 Fortran 源码逐一比对。扣除普通注释和续行格式，并将新标识符映射回旧名后，源程序的词法序列一致。导入来源检查的 250 份文件哈希核对通过。
- 按新 README 生成 V3 输入后，串行 `gj` 和双进程 `gj_mpi` 均生成与 `data/test.gj` 逐字节一致的文件，包含 7 行能级数据。
- 25 项相关 CTest 中 21 项通过，覆盖 RCI 串行和 MPI 的 Slater 积分计数、MPI Breit、串行和 MPI 的 V3 计算流程、`rwfnestimate` 头文件、串行双正交流程、Davidson、FICAL 和四种 RMCDHF 的权重。
- `state_workflow_serial`、`state_workflow_mpi`、`descriptor_input_rmcdhf_mem` 和 `descriptor_input_rmcdhf_mem_mpi` 在 macOS 上失败。两种内存版的内存探测直接读取 Linux 的 `/proc/meminfo`，本机缺少该文件。用修复前的二进制重跑，4 项均复现同一错误。全流程的内存版计算验证需要提供该文件的 Linux 环境。
- 补跑的 6 项 MCP 读取和内存缓存 CTest 全部通过，其中内存查询使用测试桩。双进程双正交和跃迁夹具也通过。该夹具在较长目录下触发了 `sys_getwd` 错误，改用工作区 `temp/tests/btfix.EXp0bj/` 后运行成功。
- 核对了 35 份 README 和本报告的 1,303 个仓库内链接、20 份应用 README 的逐文件清单，以及新示例的 fish 语法。

审查基准日志在工作区 `temp/tests/suite-source-review.gcXr3n6H/`。本轮日志在 `temp/tests/suite-source-review-fix.mXMDNWkF/`，包括 `build.log`、`ctest.log`、`ctest-mcp.log`、`fortran-comparison.log`、`gj-serial.log`、`gj-mpi.log`、`descriptor-pipeline-mpi-short-path.log` 和四份 `before-*.log`。
