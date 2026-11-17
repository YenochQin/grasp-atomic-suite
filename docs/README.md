# 文档索引与实施状态审计

最初审计日期：2026-10-05，基准 `8b298a6`；2026-10-13 更新参数恢复及严格内存 MPI 条目。
核查覆盖原 `docs/` 的全部 10 份文档，并对照当前源码、CMake、测试文件及相关
Git 提交；本次没有重跑性能基准或完整原子物理计算。

“已实施”表示对应代码/接口已经存在，不表示所有体系的精度和性能已经验收。
历史报告中的已实施功能可以继续有效，但其旧默认值、诊断输出和性能数字可能
已不适合指导当前操作。每份文档开头标明状态和适用范围。

## 目录与阅读顺序

```text
docs/
├── README.md       # 全部文档的状态、核查依据和待补事项
├── implemented/    # 已实施功能的现行说明与代码流程
├── plans/          # 待执行的代码实施与实验/验证方案
├── reviews/        # 已完成的审查与证据记录
├── reference/      # 理论背景和参考草稿
└── archive/        # 已被替代的蓝图、迁移快照和阶段性报告
```

当前运行与构建先看[仓库 README](../README_ZH.md)和
[脚本说明](../scripts/README.md)；参数管理看 `implemented/`，
代码实施计划和网格实验看 `plans/`，追溯变更和旧测量看 `archive/`。
归档文档保留设计与实验历史，不能直接当成当前配置指南。

普通串行计算链可按 [rangular 流程与共享状态](implemented/rangular_flow_and_state.md)
→ [rmcdhf 流程与共享状态](implemented/rmcdhf_flow_and_state.md) 阅读，分别对应
角系数生成及径向轨道/混合系数的 SCF 迭代；两份文档均区分现行行为与待验证的 OpenMP 候选。

## 三条独立后续任务

以下文档分别作为新对话的任务入口，包含既有结果、待办、归属、验证要求和可复制的启动提示。
本次拆分只整理文档，未启动代码实施；三个任务均默认不提交、不推送。

| 任务 | 入口文档 | 起步范围 |
| --- | --- | --- |
| 变量语义重命名 | [variable_renaming_followup_2026-11-16.md](plans/variable_renaming_followup_2026-11-16.md) | 第一批 529 项已处理，原有 19,579 项待重新对齐核实；优先保留数值过程，避免大改将退役的 CSF reader |
| CSF 描述符读取迁移 | [csf_descriptor_migration_2026-11-16.md](plans/csf_descriptor_migration_2026-11-16.md) | 2026-11-17 生产迁移与功能回归已完成；Rust 原生 Parquet/header + C ABI，全配置 CTest 通过，历史告警门禁既有失配单列 |
| 剩余非风格诊断修复 | [static_diagnostics_fixes_2026-11-16.md](plans/static_diagnostics_fixes_2026-11-16.md) | 优先 RIS4 缓存和 JJ2LSJ 唯一标签，精度/接口独立分批；CSF P1 由迁移任务承接，纯风格不动 |

各对话开始时重新核对工作树。交接时 HEAD 为 `6f96222`，已包含另一任务的 rangular OpenMP 移除；历史目标、测试数和行号不代表当前状态。
同文件的改名、输入迁移和缺陷修复须按批次排程或在独立检出中整合，不能互相覆盖。

## 逐份文档结论

| 分类 | 文档 | 实施状态与适用范围 |
| --- | --- | --- |
| 原生 CSF 输入 | [native_csf_descriptors.md](implemented/native_csf_descriptors.md) | 全部保留消费者直接读 Parquet/header；无 Python 转换和旧文本入口，含构建方式、MPI/双基组/低内存契约及验收记录。 |
| USE ONLY 第十批：PACK 接口 | [use_only_pack_2026-11-17.md](reviews/use_only_pack_2026-11-17.md) | 11 份源码明确导入重导出的 pack 过程，保留实现绑定；剩余 C121 为 1,746。 |
| USE ONLY 第十一批：RCI 与 RIS4 状态模块 | [use_only_state5_2026-11-17.md](reviews/use_only_state5_2026-11-17.md) | 29 份源码补充 33 处 ONLY，保留全部公开符号；剩余 C121 为 1,713。 |
| USE ONLY 第十二批：RIS 状态模块 | [use_only_ris_2026-11-17.md](reviews/use_only_ris_2026-11-17.md) | 14 份源码补充 15 处 ONLY；剩余 C121 为 1,698。 |
| USE ONLY 第十三批：BILST 状态模块 | [use_only_bilst_2026-11-17.md](reviews/use_only_bilst_2026-11-17.md) | 26 份源码补充 26 处 ONLY；剩余 C121 为 1,672。 |
| USE ONLY 第十四批：POS 状态模块 | [use_only_pos_2026-11-17.md](reviews/use_only_pos_2026-11-17.md) | 28 份源码补充 28 处 ONLY；剩余 C121 为 1,644。 |
| USE ONLY 第十五批：角系数缓冲区 | [use_only_buffer_2026-11-17.md](reviews/use_only_buffer_2026-11-17.md) | 15 份源码补充 15 处 ONLY；剩余 C121 为 1,629。 |
| USE ONLY 第十六批：跃迁振子强度状态模块 | [use_only_osc_2026-11-17.md](reviews/use_only_osc_2026-11-17.md) | 21 份源码补充 21 处 ONLY，剩余 C121 为 1,608。 |
| USE ONLY 第十七批：跃迁程序的 JJ2LSJ 数据导入 | [use_only_jj2lsjbio_2026-11-17.md](reviews/use_only_jj2lsjbio_2026-11-17.md) | 15 份源码补充 15 处 ONLY，剩余 C121 为 1,593。 |
| USE ONLY 第十八批：轨道变换、分块和本征解模块 | [use_only_transition_state_2026-11-17.md](reviews/use_only_transition_state_2026-11-17.md) | 138 份源码补充 164 处 ONLY，剩余 C121 为 1,429。 |
| USE ONLY 第十九批：径向网格和波函数模块 | [use_only_radial_2026-11-17.md](reviews/use_only_radial_2026-11-17.md) | 106 份源码补充 152 处 ONLY，剩余 C121 为 1,277。 |
| USE ONLY 第二十批：RMCDHF 和模块内的径向数据导入 | [use_only_radial2_2026-11-17.md](reviews/use_only_radial2_2026-11-17.md) | 107 份文件补充 163 处 ONLY，剩余 C121 为 1,114。 |
| USE ONLY 第二十一批：SCF、哈密顿矩阵和过程接口 | [use_only_scf_2026-11-17.md](reviews/use_only_scf_2026-11-17.md) | 116 份源码处理 153 处导入，剩余 C121 为 961。 |
| USE ONLY 第二十二批：分配模块与测试替身 | [use_only_memory_2026-11-17.md](reviews/use_only_memory_2026-11-17.md) | 182 份源码补充 182 处 ONLY，剩余 C121 为 779。 |
| USE ONLY 第二十三批：物理参数和轨道状态 | [use_only_def_orb_2026-11-17.md](reviews/use_only_def_orb_2026-11-17.md) | 149 份源码补充 181 处 ONLY，剩余 C121 为 598。 |
| USE ONLY 第二十四批：RMCDHF 物理参数和轨道状态 | [use_only_rmcdhf_def_orb_2026-11-17.md](reviews/use_only_rmcdhf_def_orb_2026-11-17.md) | 136 份源码补充 180 处 ONLY，剩余 C121 为 418。 |
| USE ONLY 第二十五批：轨道状态和角动量标签 | [use_only_orb_label_2026-11-17.md](reviews/use_only_orb_label_2026-11-17.md) | 174 份源码处理 180 处导入，剩余 C121 为 238。 |
| USE ONLY 第二十六批：MPI 和串行模块导入 | [use_only_mpi_2026-11-17.md](reviews/use_only_mpi_2026-11-17.md) | 143 份文件处理 166 处导入，剩余 C121 为 72。 |
| USE ONLY 第二十七批及保留例外 | [use_only_final_2026-11-17.md](reviews/use_only_final_2026-11-17.md) | 31 份文件补充 71 处 ONLY，剩余 1 条为 mpi_c 的公开 MPI 重导出。 |
| USE ONLY 第九批：数组状态 | [use_only_state_arrays_2026-11-17.md](reviews/use_only_state_arrays_2026-11-17.md) | 146 份文件补齐 193 条 ONLY，保留八个模块的完整公开实体；剩余 C121 为 1,757。 |
| USE ONLY 第八批：应用选项 | [use_only_application_flags_2026-11-17.md](reviews/use_only_application_flags_2026-11-17.md) | 138 份应用源码补齐 151 条调试/计算选项导入，剩余 C121 为 1,950。 |
| USE ONLY 第七批：调试状态 | [use_only_debug_2026-11-17.md](reviews/use_only_debug_2026-11-17.md) | 15 份共享库源码补齐 debug_c 完整导入列表，剩余 C121 为 2,101。 |
| 当前静态诊断接续入口 | [static_diagnostics_remaining_2026-11-17.md](plans/static_diagnostics_remaining_2026-11-17.md) | 提交前复查默认 6,172 条、扩展 13,604 条；USE ONLY 剩余 2,116 条（637 份文件），附位置/哈希清单，区分规范、CSF 风险和工具限制。 |
| USE ONLY 第一批：接口模块 | [use_only_interfaces_2026-11-17.md](reviews/use_only_interfaces_2026-11-17.md) | 普通四种 RMCDHF 的 142 份文件补齐 524 条 ONLY，编译器导出表确认单过程接口范围；剩余 4,365 条按依赖与重导出继续分批。 |
| USE ONLY 第二批：其他应用接口 | [use_only_app_interfaces_2026-11-17.md](reviews/use_only_app_interfaces_2026-11-17.md) | 其他 15 个应用目录的 275 份文件补齐 906 条 ONLY，247 个模块名的导出过程核对一致；剩余 3,459 条。 |
| USE ONLY 第三批：共享库过程接口 | [use_only_library_interfaces_2026-11-17.md](reviews/use_only_library_interfaces_2026-11-17.md) | 5 个共享库的 115 份文件补齐 396 条 ONLY，核对 133 个模块名，保留模块头部重导出；剩余 3,063 条。 |
| USE ONLY 第四批：小型状态模块 | [use_only_small_state_2026-11-17.md](reviews/use_only_small_state_2026-11-17.md) | 220 份文件补齐 284 条 ONLY，保留 13 个状态模块全部公开变量；剩余 2,779 条。 |
| USE ONLY 第五批：状态模块与别名 | [use_only_state_aliases_2026-11-17.md](reviews/use_only_state_aliases_2026-11-17.md) | 216 份文件补齐 383 条 ONLY，保留 15 个模块原有状态、精度别名和容量常量；剩余 2,396 条。 |
| USE ONLY 第六批：状态表与工作数组 | [use_only_state_tables_2026-11-17.md](reviews/use_only_state_tables_2026-11-17.md) | 188 份文件补齐 280 条 ONLY，保留 35 个模块的完整公开数据实体；剩余 2,116 条。 |
| 输出格式内联，待实施 | [inline_output_formats_2026-11-17.md](plans/inline_output_formats_2026-11-17.md) | 将数字标签 FORMAT 输出全部改为使用处显式格式字符串；当前统计 909 处 WRITE、781 个 FORMAT 定义，要求保持数值精度及输出字节不变。本轮仅记录。 |
| 风格第一批：空格与行尾 | [style_spacing_2026-11-17.md](reviews/style_spacing_2026-11-17.md) | 623 份文件消除 4,514 条关键字空格/行尾提示，708 个 Release 对象字节一致；有解析错误的文件跳过。 |
| 风格第二批：引号与 END 命名 | [style_quotes_end_2026-11-17.md](reviews/style_quotes_end_2026-11-17.md) | 1,287 份文件消除 14,980 条提示，字符串内容保持不变；解析问题及未提供安全修复的提示保留。 |
| 风格第三批：小数写法 | [style_decimal_2026-11-17.md](reviews/style_decimal_2026-11-17.md) | 153 份文件消除 656 条提示，保留指数和 kind；49 对常量的编译后位模式一致，解析问题文件中的 16 条提示保留。 |
| 风格第四批：解析问题文件 | [style_parser_files_2026-11-17.md](reviews/style_parser_files_2026-11-17.md) | 剩余引号、END 命名、小数提示归零；检查完整文件及 27 个预处理上下文，原始扫描的 20 条解析提示仍保留。 |
| 风格第五批：运算符、数组与空格 | [style_operators_arrays_2026-11-17.md](reviews/style_operators_arrays_2026-11-17.md) | 207 份文件统一 1,259 处关系运算符、55 处数组构造及 186 处关键字空格，相关规则归零；精度及表达式顺序不变。 |
| S1–S4 后工具全量复查 | [fortls_fortitude_rescan_2026-11-17.md](reviews/fortls_fortitude_rescan_2026-11-17.md) | 当前 1,731 份文本；Fortitude 默认 31,061 条，C181 剩 10 条；fortls 100 个上下文最终 50 条、3 条 Error 已复现为解析误报。与同工具扫描 HEAD 比较，未发现本批新增已确认功能缺陷；CSF 输入风险仍归描述符迁移任务。 |
| S3/S4 已实施，平台限制明确 | [static_numeric_interfaces_2026-11-17.md](reviews/static_numeric_interfaces_2026-11-17.md) | 20 处常量改用 real64；BNDINV 接口边界、JJ2LSJ 指针默认值和 MPI 输出失败诊断。Release/Check 软件回归通过，1,956 个加强类型检查上下文通过；普通 SCF 6 轮结果一致，Ni/Zn 初始轨道有末位差异，memory 完整 SCF 对照受 macOS 平台限制。 |
| S1/S2 首批已实施，过程级回归 | [static_input_fixes_2026-11-17.md](reviews/static_input_fixes_2026-11-17.md) | RIS4 不完整缓存安全降级与 MPI root 失败传播；JJ2LSJ 完整唯一标签过程拒绝损坏输入，正常输出与基线逐字节一致。S3/S4 见独立记录，完整 RIS4 物理结果对照未新增。 |
| 共享库过程已模块化 | [library_module_migration_2026-11-16.md](reviews/library_module_migration_2026-11-16.md) | 188 份独立接口中的 200 个过程迁入模块，保留公开 API；记录原接口属性差异、MPI 依赖拆分、测试替身和串行/MPI/OpenMP 行为回归。 |
| 剩余非风格项已分级，未实施 | [remaining_nonstyle_static_analysis_2026-11-16.md](reviews/remaining_nonstyle_static_analysis_2026-11-16.md) | 基准 9d6ee31；复核剩余输入、精度、接口与解析诊断。35 个受控输入场景定位 CSF 错误接受/零退出码、RIS4 截断缓存问题；1,107 个加强隐式类型检查的编译上下文区分声明顺序与漏声明。给出 P1–P3 建议，未修改代码。 |
| 三类优先问题已修复并回归 | [priority_input_safety_fixes_2026-11-16.md](reviews/priority_input_safety_fixes_2026-11-16.md) | MCP 状态即时检查、单 rank 故障 MPI_Abort、相关 CSF 记录检查、CPATH 分配保护及四个 EPS 的 DOUBLE 声明。Release/Check 最终各 77 项通过，非 MPI 专项 5 项通过；EPS 当前未用于执行判断，未启用阈值筛选。纯风格未改，历史告警预算缺口保留。 |
| 静态分析基线，优先项已跟进 | [fortls_fortitude_analysis_2026-11-16.md](reviews/fortls_fortitude_analysis_2026-11-16.md) | 基准 e47c59d；fortls 3.2.2 / Fortitude 0.10.1 覆盖 1,902 份文本。保留 35,361 条 Fortitude 和 57 条 fortls 的原始口径、误报与工具限制记录。三类优先项已按上条修复，其余规范和风格建议未批量处理。 |
| 明确候选已实施，CSFG 已移除 | [variable_rename_csfg_removal_2026-11-16.md](reviews/variable_rename_csfg_removal_2026-11-16.md) | 529 项候选中 528 项改名、pi 保留；746 份普通源码同步引用，原容量名保持兼容别名；移除 948 份 CSFG 源码及专用构建/测试。Release/Check 编译和各 67 项回归通过（Release 含两个维护检查的定向重跑），44 对输出字节一致；历史告警预算缺口保留。 |
| 已统一大小写，软件回归通过 | [non_csfg_lowercase_2026-11-16.md](reviews/non_csfg_lowercase_2026-11-16.md) | 排除 CSFG，1,880 份 Fortran 文本中修改 1,872 份；保留字符串、注释、宏和排版。GNU 15 Release/Check 编译通过，各 67 项回归通过（Release 含三个维护检查的定向重跑）；44 对附加输出逐字节一致。记录既有告警预算覆盖缺口和行尾空白，未做语义重命名。 |
| 明确候选已处理，其余待核实 | [non_csfg_variable_inventory_2026-11-16.md](plans/non_csfg_variable_inventory_2026-11-16.md) | 原始采集基准 caaf3a5；1,880 份 Fortran 文本、20,108 个实体。529 项明确候选已逐项标记（528 改名、pi 保留），其余 19,579 项待核实；原名、原声明行及哈希作为历史快照，当前名称与位置另有 JSON 映射。 |
| 已完成依赖与预收缩分析 | [rmcdhf_mcp_precontraction_2026-11-16.md](reviews/rmcdhf_mcp_precontraction_2026-11-16.md) | 普通 RMCDHF 的矩阵/势两条 MCP 消费路径；EOL 动态密度限制、全冻结积分收缩为 H_fixed、角向模板和迭代期收缩缓存。数学可行性与现行支持分开，未改代码或运行收缩/SCF 实验。 |
| 已完成源码分析与文件压缩实验 | [rangular_mpi_mcp_storage_2026-11-16.md](reviews/rangular_mpi_mcp_storage_2026-11-16.md) | 普通 MCP 的 rank/通道拆分、记录布局、临时峰值及 RMCDHF 读取约束；已有串行样本 Zstd 无损缩减 49.6%–64.5%，逐文件 SHA-256 一致。容器、分块压缩和冗余列号消除均未实施，生产 MPI/SCF 性能待测。 |
| 已实施，串行专项回归通过 | [rangular_angular_reuse.md](implemented/rangular_angular_reuse.md) | 默认关闭的有界 CSF 对缓存；保留顺序、占据/耦合及轨道共享关系，复用原始 T/V 系数，普通 MCP 写出不变；包含命中率、缓存容量与 RSS 测量，MPI 未接入；rangular OpenMP 已移除。 |
| 已完成输出与性能对照 | [rangular_reuse_benchmark_2026-11-14.md](reviews/rangular_reuse_benchmark_2026-11-14.md) | 独立构建修改前基准，两组数据各三轮，文本/描述符 MCP 零容差及字节一致；合成样本命中 99.93%、整程 1.21×，现有 14,551-CSF 样本几乎不命中、慢约 37%；默认保持关闭。 |
| 已实施，GNU 15 专项回归通过 | [rwfnestimate_header_input.md](implemented/rwfnestimate_header_input.md) | 普通 rwfnestimate 原生读取 rCSFs header.toml 的轨道表和含 core 总电子数；保留文本入口和 header 字节；Ni/Zn 及 TOML Release/Check 专项各 3/3，原生描述符配对哈希不变。 |
| 根因已复现并修复，物理优化待验证 | [rwfnestimate_nan_header_2026-11-12.md](reviews/rwfnestimate_nan_header_2026-11-12.md) | Ni/Zn 首次 h/i 的搜索上界变负，导致 SETPOT/TAIL 警告及 NaN 误报成功；记录数值链、suite 迁移、上游恢复；同步命名迁移及退役更新后，Release 全量 66/66 通过。 |
| 已移除被否决变体 | [rejected_rmcdhf_variants.md](implemented/rejected_rmcdhf_variants.md) | 删除 strict/diag 两个独立求解器、专用测试、构建/安装入口及当前二进制；普通四种 RMCDHF 和共享库保持原流程。 |
| 已统一普通程序命名 | [program_naming.md](implemented/program_naming.md) | 13 个保留的普通应用目录去掉 90；GJ 可执行文件改为 gj/gj_mpi，CMake、网格工具、验收和源码链接同步；重命名阶段的构建与验证记录见说明；CSFG、共享库及上游来源路径保留。 |
| 现行构建策略 | [build_profiles.md](implemented/build_profiles.md) | CSFG 已移除；生产默认 -O2，OpenMP 仅保留 jj2lsj；标准/接口/运行时诊断放 Check/Audit；serial、MPI 和 jj2lsj OpenMP 组成 24 个有效配置，16 个诊断配置沿用历史关闭 angular OpenMP 的预算 ID；历史预算范围不代表当前告警验收。 |
| 已移除，历史归档 | [rangular_openmp.md](archive/rangular_openmp.md) | 2026-11-16 基于 9d6ee31 移除 rangular_omp、并行接口、专用线程状态、构建及测试入口；串行/MPI、描述符输入和串行系数复用保留，jj2lsj OpenMP 保留。历史测量不代表当前支持范围。 |
| 已完成单机性能与输出对照 | [rangular_mpi_openmp_benchmark_2026-11-11.md](reviews/rangular_mpi_openmp_benchmark_2026-11-11.md) | 基于 c5f36d2，M4 上 14,551 CSFs、33 次 Release 运行；8 并发 MPI 4.06 秒、OpenMP 10.56 秒；描述符主要加快输入，整程收益较小；串行/线程字节一致，MPI 零容差等价。 |
| 已完成 x86 性能诊断 | [rangular_mpi_openmp_x86_diagnosis_2026-11-11.md](reviews/rangular_mpi_openmp_x86_diagnosis_2026-11-11.md) | 基于 c47f4d4，23,459 CSFs 三轮：OMP 8 线程中位数 27.06 秒、MPI 8 进程 13.49 秒；38,461,197 项贡献零差异。隔离阶段探针确认并行范围、逐列批次等待和串行结果处理的限制；未修改生产计算代码。 |
| 已实施 MPI 输入适配 | [rangular_mpi_descriptor_input.md](implemented/rangular_mpi_descriptor_input.md) | 普通 rangular_mpi 支持 rank 0 二进制加载、元数据/数组广播及逐 rank 阶段计时；保留文本入口和原 MCP 分片。Release/Check 六项 MPI 专项回归通过，提供交替重复对照脚本；不涉及 CSFG。 |
| 静态分析完成，待实施 | [mpi_csf_descriptor_input_2026-11-16.md](plans/mpi_csf_descriptor_input_2026-11-16.md) | 基准 caaf3a5；确认 rangular_mpi 已支持描述符，规划 rmcdhf_mpi/rmcdhf_mem_mpi 的占据数加载、MCP 契约检查、分块广播、错误协议及验证；保留低内存布局，单列遗留 EAL 权重依赖。未修改计算代码或运行新实验。 |
| 已完成 MPI 描述符测量 | [rangular_mpi_descriptor_benchmark_2026-11-11.md](reviews/rangular_mpi_descriptor_benchmark_2026-11-11.md) | 同一 8-rank 程序、23,459 CSFs 复核三轮：输入含广播中位数 0.411→0.061 秒，整程 12.565→12.162 秒，约减少 3.2%；96 个 MCP 文件与原 MPI 输出一致。 |
| 已实施，运行器回归通过 | [rangular_version_comparison.md](implemented/rangular_version_comparison.md) | 指定原版/新版可执行文件，独立目录运行文本及可选描述符路径；保存日志、墙钟、阶段计时与 MCP 对照；7 项维护测试通过，支持差异/超时/提前退出检测。 |
| 已实施，输入/输出回归通过 | [rangular_descriptor_input.md](implemented/rangular_descriptor_input.md) | suite 预处理二进制、Fortran reader、文本数组对照和串行阶段计时；14 项专项测试、原生 rCSFs 171 CSF 对照及小型 MCP 一致性通过；另有 14,551-CSF 串行/OpenMP 输入计时，完整科学结果及其他规模性能待验证。 |
| 首版输入和 OpenMP 已实施 | [rangular_descriptor_openmp_2026-11-10.md](plans/rangular_descriptor_openmp_2026-11-10.md) | 保留设计历史；输入、线程私有工作区、动态任务与有序提交已有首版，显式 context 和排序并行仍待后续工作。 |
| 已实施，解析器回归通过 | [mcp_compare.md](implemented/mcp_compare.md) | 2026-11-10，普通 rangular 串行/MPI MCP 解析、结构校验、按 CSF 对比较和 JSON 导出；合成回归与现有文件读回通过，不代表候选 rangular 或物理计算验收。 |
| 现行流程，静态核查完成 | [rangular_flow_and_state.md](implemented/rangular_flow_and_state.md) | 2026-11-10，基准 d181bf4；普通串行 rangular 的 block/CSF 对遍历、角系数缓冲与排序、共享工作区和文件状态。未实施 OpenMP，未运行性能或物理验证。 |
| 现行流程，静态核查完成 | [rmcdhf_flow_and_state.md](implemented/rmcdhf_flow_and_state.md) | 2026-11-10，基准 d181bf4；普通串行 rmcdhf 的矩阵/混合系数、径向积分、轨道改善、停止条件及共享状态；标出 EAL 定义路径待复核。未实施 OpenMP，未运行性能或物理验证。 |
| 已完成源码对比，数值影响待验证 | [rmcdhf_mem_comparison_2026-11-11.md](reviews/rmcdhf_mem_comparison_2026-11-11.md) | 基准 db39505；四种普通 RMCDHF 的 MCP 缓存、乘子轨道对筛选、MPI 乘子公式和径向边界差异。不能视为纯存储切换；未运行原子计算或性能基准。 |
| 外部采样首版已实现，服务器验证待进行 | [rmcdhf_mpi_diagnostics.md](implemented/rmcdhf_mpi_diagnostics.md) | 2026-11-15，普通 MPI 原版流程与兼容修复核查；Linux PID 树、D/R 线程、I/O、缺页、PSI 和 cgroup 只读采样。模拟回归通过；真实开销及停顿根因未验证。 |
| 诊断计划，D0 首版已实现 | [rmcdhf_mpi_single_node_diagnostics_2026-11-09.md](plans/rmcdhf_mpi_single_node_diagnostics_2026-11-09.md) | 普通/memory MPI 在大规模 CSF 下全程反复 D 状态；记录单节点条件、SETLAG 冻结侧冗余和 Davidson 候选。服务器采样、内部阶段计时及性能优化待进行。 |
| 优先警告修复第一批通过 | [f2018_warning_repairs_2026-11-09.md](reviews/f2018_warning_repairs_2026-11-09.md) | GETALDWT、RCI 计数、Davidson、DACON、INTERP 及不可达步长分支；正确性复核 515→475，INTENT 9→0；32 配置、1,960 次 CTest 全通过，全开 93/93。 |
| 首次矩阵/基线设施证据 | [f2018_matrix_warning_baseline_2026-11-09.md](reviews/f2018_matrix_warning_baseline_2026-11-09.md) | 保留初次 32 配置、1,560 次 CTest 通过与 3,040 文件/515 条正确性诊断记录；当前基线及修复见上一项。 |
| 已更新，剩余正确性复核待完成 | [gfortran15-warning-baseline.md](reference/f2018/gfortran15-warning-baseline.md) | 普通 Check/Audit 的 24 个配置、1,874 文件、2,200 条诊断身份，其中 337 条正确性复核身份；包含两条既有 QED 警告在本机的不同定位，不代表新增或已解决的缺陷数。 |
| Linux 初始修复证据 | [f2018_linux_davidson_restart_2026-11-09.md](reviews/f2018_linux_davidson_restart_2026-11-09.md) | 保留 Davidson、空 rank、memory IOS 修复及 71/73、strict 旧平台重现证据；strict/diag 已移除，相关测试属于历史证据。 |
| 此前严格接口验证 | [f2018_strict_interfaces_npol_2026-11-09.md](reviews/f2018_strict_interfaces_npol_2026-11-09.md) | 内部接口/NP 路径/SOLVE 意图修复；保留 macOS Check/Release 各 64/68、4 skip 及一轮全输出一致记录；后续 Linux 状态见上一项。 |
| 设施初始验证记录 | [f2018_acceptance_facilities_2026-11-09.md](reviews/f2018_acceptance_facilities_2026-11-09.md) | profiles、强制清单、模块/安装隔离与 DDOT 14/14；保留初始 audit 50/63 及严格编译失败证据，后续状态见上一项。 |
| 已移除参数放行 | [argument_interfaces_2026-11-08.md](reviews/argument_interfaces_2026-11-08.md) | gfortran 15/F2018 全量构建不再使用 fallow-argument-mismatch；MPI/BLAS 显式接口、类型和有界数组区段修复；当前 CTest 36/47，11 项为既有失败；orbopt 已退役。 |
| 已移除旧程序 | [orbopt_retirement.md](implemented/orbopt_retirement.md) | 2026-11-08 删除四个 rmcdhf_orbopt 目标/源码和专用驱动，记录 2026-11-08 的 orbopt 退役；strict/diag 后续也已移除，现行网格入口为 19 个。 |
| 已实施兼容性子批次 | [f2018_compatibility.md](implemented/f2018_compatibility.md) | 首轮标准语法及后续显式接口、NP 路径修复；严格 Check/Release 全量编译/链接通过，完整数值及静态状态验收未完成。 |
| 已完成本轮验证 | [f2018_gfortran15_batches_2026-11-08.md](reviews/f2018_gfortran15_batches_2026-11-08.md) | 全量构建通过；五项新增回归通过；四个 solver 一轮 SCF 共 28 文件与旧版本一致；CTest 33/44，11 个失败与旧基线相同，无新增失败。 |
| 已完成首轮编译诊断 | [f2018_gfortran15_compile_probe_2026-11-08.md](reviews/f2018_gfortran15_compile_probe_2026-11-08.md) | GNU 构建已加 `-std=f2018`，全新全量构建失败；15 个文件直接报错。补充生产源码诊断为 59 个文件失败，移除参数放行后为 143 个；列出共享库、平台函数、kind/语法与 MPI/BLAS 接口优先项，未修改计算实现或运行数值测试。 |
| 已实施 | [standard_grasp_integration.md](implemented/standard_grasp_integration.md) | 2026-11-08 导入普通 rangular、rwfnestimate、rbiotransform、rci、rtransition 五类程序共十个目标，复用 suite 共享库/参数；移除初猜程序局部容量定义，补丁器覆盖 25 个网格/4 个续算入口。完整 F2018 与原子物理验收尚未完成。 |
| 已完成接入检查 | [standard_grasp_integration_validation_2026-11-08.md](reviews/standard_grasp_integration_validation_2026-11-08.md) | GNU 15.3.0 全量及 serial-only 构建通过；35 项补丁器测试通过；CTest 28/39 通过，11 项 NP/strict 失败在未修改基线复现。记录既有并行模块输出冲突；不宣称全绿或物理验收。 |
| 已实施 | [jj2lsj_postprocessing.md](implemented/jj2lsj_postprocessing.md) | 2026-11-04 分阶段计时、稳定归并排序及 ASF 排序并行；唯一标签和输出按序执行，跨 block 流水线经计时评估暂缓。 |
| 已实施 | [jj2lsj_openmp.md](implemented/jj2lsj_openmp.md) | 2026-11-02 `jj2lsj_omp` 普通 CSF 的 OpenMP LS 行转换，保留默认提前停止与串行 `.lsj.T`；数值/文件回归及 Debug、无 OpenMP 构建通过；Gd I 真实输入兼容性通过，已有单次性能观测。 |
| 已移除，历史归档 | [graspg_integration.md](archive/graspg_integration.md) | 记录 2026-10-30 的 GRASPG 集成与来源；2026-11-16 专用源码、构建和测试入口已移除。旧启用说明不是当前可用功能，普通共享库保留。 |
| 已实施 | [rdensity_mpi.md](implemented/rdensity_mpi.md) | 2026-10-09 新增串行/MPI 自然轨道程序；CSF 分工、矩阵归约、解析小模型及 Ni 对照测试。大规模性能与物理收敛仍待验证。 |
| 已实施 | [common_parameters.md](implemented/common_parameters.md) | 2026-11-16 使用语义容量名称并保留旧名别名；CSFG 移除后覆盖 15 个普通网格入口、2 个 RCI 续算入口，默认容量仍为 590/600。 |
| 已实施 | [RHFS_gJ_report.md](implemented/RHFS_gJ_report.md) | 串行 RHFS 的加载、矩阵元、ASF 投影及输出链仍存在。仅覆盖该路径，不是 gj/MPI 完整说明。 |
| 已实施，普通 MPI 已移除 NP 入口 | [rmcdhf_nuclear_polarization.md](implemented/rmcdhf_nuclear_polarization.md) | 2026-11-15 普通 rmcdhf_mpi 恢复直接核势初始化；rmcdhf、rmcdhf_mem、rmcdhf_mem_mpi 保留可选 NP。保留 2026-10-24 历史代码验证，物理计算验证未执行。 |
| 已完成代码计划 | [rmcdhf_nuclear_polarization_implementation.md](plans/rmcdhf_nuclear_polarization_implementation.md) | 保留实施设计、完成清单与结果；现行输入契约见已实施说明。 |
| 待执行 | [grasp_grid_parameter_tuning.md](plans/grasp_grid_parameter_tuning.md) | 修改脚本已实现，文档配置已核查；D/L/O/ACCY 和生产模型物理验收仍待执行。 |
| 优先修复已推进，完整迁移未完成 | [f2018_migration_gfortran15_2026-11-08.md](plans/f2018_migration_gfortran15_2026-11-08.md) | 当前 32 配置、1,960 次 CTest 通过；正确性复核诊断减少 40，基线显式收紧。剩余警告、状态/接口、真实多根与物理验收待完成。 |
| 阶段性实施 | [jj2lsj_postprocessing_parallelization_2026-11-03.md](plans/jj2lsj_postprocessing_parallelization_2026-11-03.md) | 后处理按 ASF 并行、唯一标签依赖与跨 block 异步流水线讨论；计时与稳定排序已实施，保留 ASF 后处理和跨 block 流水线的后续设计。 |
| 已完成验证 | [jj2lsj_postprocessing_validation_2026-11-04.md](reviews/jj2lsj_postprocessing_validation_2026-11-04.md) | AS4 稳定排序及整理约 9.60→0.045 s，总耗时单次约减少 2.1%；AS4/AS1 输出一致，标签异步化经计时评估暂缓。 |
| 已完成验证 | [jj2lsj_gd_i_validation_2026-11-02.md](reviews/jj2lsj_gd_i_validation_2026-11-02.md) | 六套 Gd I 输入的 18 个输出文件一致；两个代表案例的 994,894 个 LS 系数在 4/8 线程下与串行逐位一致，包含真实提前停止路径。 |
| 已完成审查 | [grasp_grid_methodology_review.md](reviews/grasp_grid_methodology_review.md) | 源码审查和最小例程诊断已完成，文档建议已落实。旧 P 配置等问题描述有指定历史版本。 |
| 理论参考 | [lande_g.md](reference/lande_g.md) | 概念、ASF 与 LS 近似的参考草稿，不能划为功能实施计划。算符规范与系数尚未完整独立核对。 |
| 参数参考 | [nuclear_polarization_parameters.md](reference/nuclear_polarization_parameters.md) | 2026-10-24 核对 Flambaum 2021 的 E1/E2 强度、50 keV 修正极化率、截断参数与原子单位换算；记录拟合系数内部差异，未建立生产核素参数库或物理精度结论。 |
| 历史蓝图 | [gj90_blueprint.md](archive/gj90_blueprint.md) | 阶段 A/C 主体已实现，B 有历史样例对比记录；阶段 D 的 OpenMP/线程重构未见实现。不再是当前开发清单。 |
| 历史阶段报告 | [RHFS_MPI_hfsgg_optimization.md](archive/RHFS_MPI_hfsgg_optimization.md) | 前置筛选仍在；诊断计数、直接 ICHKQ1 调用和“当前瓶颈”描述已被后续代码覆盖。 |
| 历史分支总结 | [1.0.1-dev.1_optimization_summary.md](archive/1.0.1-dev.1_optimization_summary.md) | HFSGG 与 DIAGA3/4 保留。RECOP2 独立快速预检查未保留，已修正原文误记；旧计时未复测。 |
| 历史分支总结 | [1.1.1-dev.1_optimization_summary.md](archive/1.1.1-dev.1_optimization_summary.md) | RIS 主要优化仍存在；历史耗时、测试名称和测试覆盖不代表当前数值验收。 |
| 历史迁移记录 | [rmcdhf_migration.md](archive/rmcdhf_migration.md) | 四个 orbopt 目标和共享库已落地；容量经转导层公开，默认值和网格初始化恢复原版，以现行参数说明为准。 |

## 核查依据与具体过时项

### 参数与仓库迁移

- [suite_parameters_M.f90](../src/lib/libmod/suite_parameters_M.f90) 定义容量，
  [parameter_def_M.f90](../src/lib/libmod/parameter_def_M.f90) 转导原有常量名；
  当前 `NNNP=590`、`NNN1=600`。网格初始化仍在各入口。
- 当前网格脚本覆盖十五个普通网格入口和两个 RCI 续算 ACCY 入口，保留各自的默认值与初始化时序；
  CSFG 和此前退役的求解器不在现行范围。有限核默认 `radial_point_count=max_radial_points=590`，旧容量名仍由转导接口提供。
- [test/CMakeLists.txt](../test/CMakeLists.txt) 保留共享容量检查和网格脚本回归；
  历史测试数量不能代替当前实际验证；退役求解器的测试记录仅供追溯。
- 历史 `N=2990,H=0.05` 的风险审查与候选网格实验仍有参考价值，
  但该组合已不是仓库默认值，参数恢复本身也不代表生产体系收敛。

### RHFS 与独立 gj

- [hfs92.f90](../src/appl/rhfs/hfs92.f90) 与
  [hfsgg.f90](../src/appl/rhfs/hfsgg.f90) 对应现行串行流程报告。
- [gj.f90](../src/appl/gj/gj.f90)、
  [compute_gj_only.f90](../src/appl/gj/compute_gj_only.f90) 和
  [setout_gj.f90](../src/appl/gj/setout_gj.f90) 已实现独立入口、NVEC 对角累加
  和 `.gj/.cgj` 输出，取代蓝图的“待新建”描述。
- [compute_gj_only_mpi.f90](../src/appl/gj/compute_gj_only_mpi.f90) 已按
  `IC=MYID+1,NCF,NPROCS` 分工并归约；主循环 OpenMP/线程工作区重构未见实现。
- 蓝图曾把 RHFS `getmixblock.f90` 写在 `lib9290/` 下，实际文件位于
  [src/appl/rhfs](../src/appl/rhfs/getmixblock.f90)，归档正文已更正路径。

### RHFS/RIS 性能优化

- [hfsgg_mpi.f90](../src/appl/rhfs/hfsgg_mpi.f90) 保留三角/宇称、内联 IQA
  占据数筛选，以及 ELEMFAC/GJFAC/DGJFAC 的预计算和零贡献筛选。
  早期 `cand/ptrig/zelem` 等 profiling 输出已在 `3d067e7` 删除。
- `f414fb6` 增加的 RECOP2 独立 `IRE=0` 分支在 `505d7f4` 中删除，当前
  [recop2.f90](../src/lib/librang90/recop2.f90) 没有该分支。
  [diaga3.f90](../src/lib/librang90/diaga3.f90) 与
  [diaga4.f90](../src/lib/librang90/diaga4.f90) 的快速预检查仍在。
- [ris_cal_mpi.f90](../src/appl/ris4/ris_cal_mpi.f90) 的 YES2 MPI 分支和
  `GDRSUMMPI_ROOT`，以及密度/SMS 内核的外层筛选与 EVPAIR 复用仍在。
  两份分支总结的性能数字只说明当时样例，不可作为当前硬件/模型的速度承诺。

## 仍需补齐的验证

1. 网格方案的 Ni/Cl/U 或所选代表体系及最终生产模型收敛计算，按方案逐阶段执行。
2. 当前 [CMakeLists.txt](../CMakeLists.txt) 的 gj 数值回归注册依赖
   `tests/gj_regression.sh` 和 `tests/gj_mpi_regression.sh`，但当前受版本控制
   的源码树没有这些脚本，因而这些条件不会注册对应测试。
   不能用旧分支报告的 CTest 记录证明当前自动 g_J 回归已覆盖；需要补齐脚本或
   调整测试接入。本次分类只记录该缺口，没有修改构建实现。
3. RIS 优化总结明确没有专门的 RIS 数值回归；需要实际生成/读回路径与串行/MPI
   对照结果，不能只以 g_J 或公共库 CTest 通过代替。
4. 理论参考稿的完整算符、单位及约化矩阵元约定仍需第一手文献核对。

保留求解器的共有检查见[测试入口](../test/rmcdhf/README.md)；旧 orbopt 驱动说明和
[历史结果](archive/orbopt/RESULTS.md)已归档。另保留
[已退役的 stage guard 历史记录](archive/orbopt_stage_guard_retired.md)，
2026-10-19 移除 Python 阶段验收流程及其专属入口，实施记录已归档；
上述记录描述当时状态；2026-11-08 又移除了四个 orbopt 求解器及其专用测试驱动。
没有迁入大型外部计算数据。

## 后续维护规则

新增文档先按用途选择这五类目录，写明日期、代码版本、已实施部分和待验证部分。
完成计划时更新状态并指向实现证据；归档时保留历史参数和测量背景，注明当前替代入口。
源码变更触及现行说明时同步更新文档与此索引。仓库内部链接使用相对路径，
避免依赖某台机器的绝对目录和已经漂移的源码行号。
