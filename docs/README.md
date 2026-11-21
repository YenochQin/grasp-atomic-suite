# 文档索引

2026-11-20 更新输出格式内联记录。
“已实施”表示代码/接口存在，不表示所有体系的精度、收敛和性能已验收。

## 阅读入口

- 构建与运行：[仓库说明](../README_ZH.md)、[构建策略](implemented/build_profiles.md)、
  [Rust 原生 CSF 输入](implemented/native_csf_descriptors.md)、[脚本说明](../scripts/README.md)。
- 容量与网格：[共同参数](implemented/common_parameters.md)；参数容量与运行时网格点数分别核定。
- 普通计算链：[rangular 流程](implemented/rangular_flow_and_state.md) →
  [rmcdhf 流程](implemented/rmcdhf_flow_and_state.md)。
- 后续代码工作：[变量语义重命名](plans/variable_renaming_followup_2026-11-16.md)、
  [静态诊断剩余事项](plans/static_diagnostics_remaining_2026-11-17.md)。

`implemented/` 保留功能契约，`plans/` 保留未完成工作，`reviews/` 保留独立分析和验证证据，
`reference/` 保存理论及机器可读来源，`archive/` 保存仍有参考价值的历史背景。
旧输入说明和退役程序记录按其标注版本阅读，不能据此恢复旧生产入口。

## 现行功能说明与保留的实现背景

| 文档 | 范围与状态 |
| --- | --- |
| [inline_output_formats.md](implemented/inline_output_formats.md) | 数字标签 `FORMAT` 内联，覆盖共享库和应用的 903 处 `WRITE`、2 处 `PRINT`；源码、预处理和字节对照已验证。 |
| [shared_library_naming.md](implemented/shared_library_naming.md) | 共享状态与基础接口库的命名约定、Fortran 调用方迁移和完整模块对照。 |
| [shared_library_procedure_names.md](implemented/shared_library_procedure_names.md) | 共享库过程名称迁移表与现行 Fortran 调用入口。 |
| [native_csf_descriptors.md](implemented/native_csf_descriptors.md) | 现行生产 CSF 输入；reader 已内置于 suite 并构建到 lib/，无需外部 rCSFs reader；所有保留消费者直接读取 Parquet/header，含 MPI、双基组、低内存及错误传播契约。 |
| [rangular_angular_reuse.md](implemented/rangular_angular_reuse.md) | 已实施，串行专项回归通过；默认关闭的有界 CSF 对缓存；保留顺序、占据/耦合及轨道共享关系，复用原始 T/V 系数，普通 MCP 写出不变；包含命中率、缓存容量与 RSS 测量，MPI 未接入；rangular OpenMP 已移除。 |
| [rwfnestimate_header_input.md](implemented/rwfnestimate_header_input.md) | header 输入说明；保留 header 轨道表、电子数和初猜检查背景；当前仅 header 输入，无旧文本回退。 |
| [rejected_rmcdhf_variants.md](implemented/rejected_rmcdhf_variants.md) | 已移除被否决变体；删除 strict/diag 两个独立求解器、专用测试、构建/安装入口及当前二进制；普通四种 RMCDHF 和共享库保持原流程。 |
| [program_naming.md](implemented/program_naming.md) | 普通程序命名；保留应用目录和可执行文件去掉 90 的命名契约；CSFG 已随后移除，共享库及上游来源命名保留。 |
| [build_profiles.md](implemented/build_profiles.md) | 现行构建策略；CSFG 已移除；生产默认 -O2，OpenMP 仅保留 jj2lsj；标准/接口/运行时诊断放 Check/Audit；serial、MPI 和 jj2lsj OpenMP 组成 24 个有效配置，16 个诊断配置沿用历史关闭 angular OpenMP 的预算 ID；历史预算范围不代表当前告警验收。 |
| [rangular_mpi_descriptor_input.md](implemented/rangular_mpi_descriptor_input.md) | 旧 MPI 输入记录；保留 rank 0 加载、广播及旧输入计时背景；GCSDESC1 入口已由 Rust 原生输入替换。 |
| [rangular_version_comparison.md](implemented/rangular_version_comparison.md) | 已实施，运行器回归通过；指定原版/新版可执行文件，独立目录运行文本及可选描述符路径；保存日志、墙钟、阶段计时与 MCP 对照；7 项维护测试通过，支持差异/超时/提前退出检测。 |
| [rangular_descriptor_input.md](implemented/rangular_descriptor_input.md) | 旧输入实现记录；旧文本/GCSDESC1 契约及数组对照方法；当前生产输入统一见 native_csf_descriptors.md。 |
| [mcp_compare.md](implemented/mcp_compare.md) | 已实施，解析器回归通过；2026-11-10，普通 rangular 串行/MPI MCP 解析、结构校验、按 CSF 对比较和 JSON 导出；合成回归与现有文件读回通过，不代表候选 rangular 或物理计算验收。 |
| [rangular_flow_and_state.md](implemented/rangular_flow_and_state.md) | 现行流程，静态核查完成；2026-11-10，基准 d181bf4；普通串行 rangular 的 block/CSF 对遍历、角系数缓冲与排序、共享工作区和文件状态。未实施 OpenMP，未运行性能或物理验证。 |
| [rmcdhf_flow_and_state.md](implemented/rmcdhf_flow_and_state.md) | 现行流程，静态核查完成；2026-11-10，基准 d181bf4；普通串行 rmcdhf 的矩阵/混合系数、径向积分、轨道改善、停止条件及共享状态；标出 EAL 定义路径待复核。未实施 OpenMP，未运行性能或物理验证。 |
| [rmcdhf_mpi_diagnostics.md](implemented/rmcdhf_mpi_diagnostics.md) | 外部采样首版已实现，服务器验证待进行；2026-11-15，普通 MPI 原版流程与兼容修复核查；Linux PID 树、D/R 线程、I/O、缺页、PSI 和 cgroup 只读采样。模拟回归通过；真实开销及停顿根因未验证。 |
| [orbopt_retirement.md](implemented/orbopt_retirement.md) | 已移除旧程序；2026-11-08 删除四个 rmcdhf_orbopt 目标/源码和专用驱动，记录 2026-11-08 的 orbopt 退役；strict/diag 后续也已移除，当前网格范围以共同参数说明为准。 |
| [f2018_compatibility.md](implemented/f2018_compatibility.md) | 现行兼容性与接口；汇总标准 kind、共享库模块 API、MPI 重导出例外、输入/数值修复、来源校验与验证限制。 |
| [core_naming.md](implemented/core_naming.md) | 公共基础库的目录、静态库、Fortran 文件与模块命名对照。 |
| [Core 计算公式](../src/lib/core/readme.md) | 公共基础库中的径向网格、求积、轨道、核势、角动量和矩阵运算公式。 |
| [core_variable_migration.md](implemented/core_variable_migration.md) | `libcore` 变量与常量的逐作用域对应表，包含形参、保留名称与共享状态别名。 |
| [breit_angular_support_migration.md](implemented/breit_angular_support_migration.md) | Breit 角系数支持库的命名迁移、共享径向积分系数过程和角系数收集接口。 |
| [spin_angular_naming.md](implemented/spin_angular_naming.md) | 自旋–角向库的目录、库目标、89 个源文件和 Fortran 模块命名对照。 |
| [librang90_variable_migration.md](implemented/librang90_variable_migration.md) | `spin_angular` 变量与常量的逐过程对照，涵盖子壳层位置、角动量编码、张量矩阵元和 Breit 角向系数。 |
| [standard_grasp_integration.md](implemented/standard_grasp_integration.md) | 普通程序接入；普通角系数、初猜、RCI 和跃迁程序复用套件共享库；历史接入检查不代表完整科学验收。 |
| [rlevels.md](implemented/rlevels.md) | 原版 `rlevels` 的 Suite 接入，含 `.m`、`.cm`、LSJ 标签匹配规则和回归入口。 |
| [jj2lsj_postprocessing.md](implemented/jj2lsj_postprocessing.md) | 已实施；2026-11-04 分阶段计时、稳定归并排序及 ASF 排序并行；唯一标签和输出按序执行，跨 block 流水线经计时评估暂缓。 |
| [jj2lsj_openmp.md](implemented/jj2lsj_openmp.md) | 已实施；2026-11-02 `jj2lsj_omp` 普通 CSF 的 OpenMP LS 行转换，保留默认提前停止与串行 `.lsj.T`；数值/文件回归及 Debug、无 OpenMP 构建通过；Gd I 真实输入兼容性通过，已有单次性能观测。 |
| [rdensity_mpi.md](implemented/rdensity_mpi.md) | 已实施；2026-10-09 新增串行/MPI 自然轨道程序；CSF 分工、矩阵归约、解析小模型及 Ni 对照测试。大规模性能与物理收敛仍待验证。 |
| [common_parameters.md](implemented/common_parameters.md) | 已实施；2026-11-16 使用语义容量名称并保留旧名别名；CSFG 移除后覆盖 15 个普通网格入口、2 个 RCI 续算入口，默认容量仍为 590/600。 |
| [davidson_eigensolver_migration.md](implemented/davidson_eigensolver_migration.md) | Davidson 求解库的目录、静态库、Fortran 文件、模块和变量命名对照。 |
| [grasp_mpi_naming.md](implemented/grasp_mpi_naming.md) | MPI 支持库的目录、静态库、Fortran 文件和模块命名对照。 |
| [mpi90_variable_migration.md](implemented/mpi90_variable_migration.md) | MPI 支持库的作用域级变量与常量对应表，含通信、输入、目录管理和稀疏矩阵例程。 |
| [RHFS_gJ_report.md](implemented/RHFS_gJ_report.md) | 已实施；串行 RHFS 的加载、矩阵元、ASF 投影及输出链仍存在。仅覆盖该路径，不是 gj/MPI 完整说明。 |
| [rmcdhf_nuclear_polarization.md](implemented/rmcdhf_nuclear_polarization.md) | 已实施，普通 MPI 已移除 NP 入口；2026-11-15 普通 rmcdhf_mpi 恢复直接核势初始化；rmcdhf、rmcdhf_mem、rmcdhf_mem_mpi 保留可选 NP。保留 2026-10-24 历史代码验证，物理计算验证未执行。 |

## 未完成的计划与清单

| 文档 | 范围与状态 |
| --- | --- |
| [CSF 计算描述符 V3](plans/csf_descriptor_v3_2026-11-21.md) | 待实施。按轨道读取完整角动量码，分别读取总 2J+1 与宇称；覆盖角积分、RCI、RMCDHF 和 MPI 消费路径。 |
| [variable_renaming_followup_2026-11-16.md](plans/variable_renaming_followup_2026-11-16.md) | 变量语义重命名；第一批 529 项已处理，原有 19,579 项需与当前源码重新对齐；旧 CSF reader 已替换，不恢复退役代码。 |
| [static_diagnostics_remaining_2026-11-17.md](plans/static_diagnostics_remaining_2026-11-17.md) | 待核定的静态诊断；S1–S4、格式及 USE ONLY 批次已完成；保留接口、状态、I/O、工具覆盖、告警预算及数值验证待办。 |
| [non_csfg_variable_inventory_2026-11-16.md](plans/non_csfg_variable_inventory_2026-11-16.md) | 明确候选已处理，其余待核实；原始采集基准 caaf3a5；1,880 份 Fortran 文本、20,108 个实体。529 项明确候选已逐项标记（528 改名、pi 保留），其余 19,579 项待核实；原名、原声明行及哈希作为历史快照，当前名称与位置另有 JSON 映射。 |
| [rmcdhf_mpi_single_node_diagnostics_2026-11-09.md](plans/rmcdhf_mpi_single_node_diagnostics_2026-11-09.md) | 诊断计划，D0 首版已实现；普通/memory MPI 在大规模 CSF 下全程反复 D 状态；记录单节点条件、SETLAG 冻结侧冗余和 Davidson 候选。服务器采样、内部阶段计时及性能优化待进行。 |
| [grasp_grid_parameter_tuning.md](plans/grasp_grid_parameter_tuning.md) | 待执行；修改脚本已实现，文档配置已核查；D/L/O/ACCY 和生产模型物理验收仍待执行。 |
| [f2018_migration_gfortran15_2026-11-08.md](plans/f2018_migration_gfortran15_2026-11-08.md) | 部分完成，验收设计保留；已完成编译和接口子批次；保留真实多根、数值、状态与平台验收设计，旧目标及测试数明确按历史版本阅读。 |
| [jj2lsj_postprocessing_parallelization_2026-11-03.md](plans/jj2lsj_postprocessing_parallelization_2026-11-03.md) | 阶段性实施；后处理按 ASF 并行、唯一标签依赖与跨 block 异步流水线讨论；计时与稳定排序已实施，保留 ASF 后处理和跨 block 流水线的后续设计。 |

## 分析与验证证据

| 文档 | 范围与状态 |
| --- | --- |
| [fortran_arrow_audit_2026-11-20.md](reviews/fortran_arrow_audit_2026-11-20.md) | Fortran 导入重命名与旧命名兼容入口已清除，含迁移对照、旧接口拒绝验证及运行结果比较。 |
| [native_csf_descriptor_validation_2026-11-17.md](reviews/native_csf_descriptor_validation_2026-11-17.md) | 原生 CSF 输入的功能回归、矩阵及告警门禁边界；保留独立验收证据。 |
| [fortls_fortitude_current_audit_2026-11-17.md](reviews/fortls_fortitude_current_audit_2026-11-17.md) | 迁移前工具审查；基准 5806884；保留原静态扫描及解析边界。之后生产 CSF reader 已替换，不能将旧风险与计数当作当前状态。 |
| [fortls_fortitude_rescan_2026-11-17.md](reviews/fortls_fortitude_rescan_2026-11-17.md) | 工具对照证据；保留 S1–S4 后的解析器/编译器对照与检查边界；旧计数和 CSF 输入风险已由后续修改覆盖，不能直接作为当前执行清单。 |
| [rmcdhf_mcp_precontraction_2026-11-16.md](reviews/rmcdhf_mcp_precontraction_2026-11-16.md) | 已完成依赖与预收缩分析；普通 RMCDHF 的矩阵/势两条 MCP 消费路径；EOL 动态密度限制、全冻结积分收缩为 H_fixed、角向模板和迭代期收缩缓存。数学可行性与现行支持分开，未改代码或运行收缩/SCF 实验。 |
| [rangular_mpi_mcp_storage_2026-11-16.md](reviews/rangular_mpi_mcp_storage_2026-11-16.md) | 已完成源码分析与文件压缩实验；普通 MCP 的 rank/通道拆分、记录布局、临时峰值及 RMCDHF 读取约束；已有串行样本 Zstd 无损缩减 49.6%–64.5%，逐文件 SHA-256 一致。容器、分块压缩和冗余列号消除均未实施，生产 MPI/SCF 性能待测。 |
| [rangular_reuse_benchmark_2026-11-14.md](reviews/rangular_reuse_benchmark_2026-11-14.md) | 已完成输出与性能对照；独立构建修改前基准，两组数据各三轮，文本/描述符 MCP 零容差及字节一致；合成样本命中 99.93%、整程 1.21×，现有 14,551-CSF 样本几乎不命中、慢约 37%；默认保持关闭。 |
| [rwfnestimate_nan_header_2026-11-12.md](reviews/rwfnestimate_nan_header_2026-11-12.md) | 根因已复现并修复，物理优化待验证；Ni/Zn 首次 h/i 的搜索上界变负，导致 SETPOT/TAIL 警告及 NaN 误报成功；记录数值链、suite 迁移、上游恢复；同步命名迁移及退役更新后，Release 全量 66/66 通过。 |
| [rangular_mpi_openmp_benchmark_2026-11-11.md](reviews/rangular_mpi_openmp_benchmark_2026-11-11.md) | 已完成单机性能与输出对照；基于 c5f36d2，M4 上 14,551 CSFs、33 次 Release 运行；8 并发 MPI 4.06 秒、OpenMP 10.56 秒；描述符主要加快输入，整程收益较小；串行/线程字节一致，MPI 零容差等价。 |
| [rangular_mpi_openmp_x86_diagnosis_2026-11-11.md](reviews/rangular_mpi_openmp_x86_diagnosis_2026-11-11.md) | 已完成 x86 性能诊断；基于 c47f4d4，23,459 CSFs 三轮：OMP 8 线程中位数 27.06 秒、MPI 8 进程 13.49 秒；38,461,197 项贡献零差异。隔离阶段探针确认并行范围、逐列批次等待和串行结果处理的限制；未修改生产计算代码。 |
| [rangular_mpi_descriptor_benchmark_2026-11-11.md](reviews/rangular_mpi_descriptor_benchmark_2026-11-11.md) | 已完成 MPI 描述符测量；同一 8-rank 程序、23,459 CSFs 复核三轮：输入含广播中位数 0.411→0.061 秒，整程 12.565→12.162 秒，约减少 3.2%；96 个 MCP 文件与原 MPI 输出一致。 |
| [rmcdhf_mem_comparison_2026-11-11.md](reviews/rmcdhf_mem_comparison_2026-11-11.md) | 已完成源码对比，数值影响待验证；基准 db39505；四种普通 RMCDHF 的 MCP 缓存、乘子轨道对筛选、MPI 乘子公式和径向边界差异。不能视为纯存储切换；未运行原子计算或性能基准。 |
| [jj2lsj_postprocessing_validation_2026-11-04.md](reviews/jj2lsj_postprocessing_validation_2026-11-04.md) | 已完成验证；AS4 稳定排序及整理约 9.60→0.045 s，总耗时单次约减少 2.1%；AS4/AS1 输出一致，标签异步化经计时评估暂缓。 |
| [jj2lsj_gd_i_validation_2026-11-02.md](reviews/jj2lsj_gd_i_validation_2026-11-02.md) | 已完成验证；六套 Gd I 输入的 18 个输出文件一致；两个代表案例的 994,894 个 LS 系数在 4/8 线程下与串行逐位一致，包含真实提前停止路径。 |
| [grasp_grid_methodology_review.md](reviews/grasp_grid_methodology_review.md) | 已完成审查；源码审查和最小例程诊断已完成，文档建议已落实。旧 P 配置等问题描述有指定历史版本。 |

## 参考资料

| 文档 | 范围与状态 |
| --- | --- |
| [gfortran15-warning-baseline.md](reference/f2018/gfortran15-warning-baseline.md) | 已更新，剩余正确性复核待完成；普通 Check/Audit 的 24 个配置、1,874 文件、2,200 条诊断身份，其中 337 条正确性复核身份；包含两条既有 QED 警告在本机的不同定位，不代表新增或已解决的缺陷数。 |
| [lande_g.md](reference/lande_g.md) | 理论参考；概念、ASF 与 LS 近似的参考草稿，不能划为功能实施计划。算符规范与系数尚未完整独立核对。 |
| [nuclear_polarization_parameters.md](reference/nuclear_polarization_parameters.md) | 参数参考；2026-10-24 核对 Flambaum 2021 的 E1/E2 强度、50 keV 修正极化率、截断参数与原子单位换算；记录拟合系数内部差异，未建立生产核素参数库或物理精度结论。 |

## 历史背景

| 文档 | 范围与状态 |
| --- | --- |
| [rangular_openmp.md](archive/rangular_openmp.md) | 已移除，历史归档；2026-11-16 基于 9d6ee31 移除 rangular_omp、并行接口、专用线程状态、构建及测试入口；串行/MPI、描述符输入和串行系数复用保留，jj2lsj OpenMP 保留。历史测量不代表当前支持范围。 |
| [graspg_integration.md](archive/graspg_integration.md) | 已移除，历史归档；记录 2026-10-30 的 GRASPG 集成与来源；2026-11-16 专用源码、构建和测试入口已移除。旧启用说明不是当前可用功能，普通共享库保留。 |
| [gj90_blueprint.md](archive/gj90_blueprint.md) | 历史蓝图；阶段 A/C 主体已实现，B 有历史样例对比记录；阶段 D 的 OpenMP/线程重构未见实现。不再是当前开发清单。 |
| [RHFS_MPI_hfsgg_optimization.md](archive/RHFS_MPI_hfsgg_optimization.md) | 历史阶段报告；前置筛选仍在；诊断计数、直接 ICHKQ1 调用和“当前瓶颈”描述已被后续代码覆盖。 |
| [1.0.1-dev.1_optimization_summary.md](archive/1.0.1-dev.1_optimization_summary.md) | 历史分支总结；HFSGG 与 DIAGA3/4 保留。RECOP2 独立快速预检查未保留，已修正原文误记；旧计时未复测。 |
| [1.1.1-dev.1_optimization_summary.md](archive/1.1.1-dev.1_optimization_summary.md) | 历史分支总结；RIS 主要优化仍存在；历史耗时、测试名称和测试覆盖不代表当前数值验收。 |
| [rmcdhf_migration.md](archive/rmcdhf_migration.md) | 历史迁移记录；四个 orbopt 目标和共享库已落地；容量经转导层公开，默认值和网格初始化恢复原版，以现行参数说明为准。 |

补充历史资料：[orbopt 驱动说明](archive/orbopt/README.md)、
[结果](archive/orbopt/RESULTS.md)、[已退役的 stage guard](archive/orbopt_stage_guard_retired.md)、
[GRASPG 上游说明](reference/graspg/UPSTREAM_README.md)。

## 验证边界与维护

网格及活性空间收敛、真实多根/完整性质计算、RIS 数值对照、服务器 MPI 性能及理论算符约定，
按各自计划与说明中的未完成范围继续核定。历史测试数不代表当前自动回归覆盖。
现有 gj 回归脚本缺口及 RIS 专项数值覆盖问题不能靠旧优化报告关闭；测试入口见
[test/CMakeLists.txt](../test/CMakeLists.txt)和[求解器测试说明](../test/rmcdhf/README.md)。

已完成的 USE ONLY、格式整理、大小写/首批改名、模块化、输入/数值修复和 F2018 分批报告，
以及完成的 CSF/NP 代码计划不再逐批保留；实现约束集中在
[兼容性说明](implemented/f2018_compatibility.md)、原生输入及 NP 说明，原始记录由 Git 历史追溯。
测试引用的来源快照、告警基线和仍供后续工作使用的清单继续保留。

后续计划完成时，将仍适用的接口、操作方式及验证限制并入现行说明，删除重复的过程记录并同步引用。
仅当报告包含独立数值证据、未解决事项或必要历史背景时继续保留；写明日期和源码基准。
本机日志、大型输入/输出和环境不纳入文档提交。
