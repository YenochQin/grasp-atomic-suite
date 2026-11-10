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

## 逐份文档结论

| 分类 | 文档 | 实施状态与适用范围 |
| --- | --- | --- |
| 现行构建策略 | [build_profiles.md](implemented/build_profiles.md) | CSFG 默认隔离；生产默认 -O2 与目标所需 OpenMP；标准/接口/运行时诊断放 Check/Audit，移除 legacy/隐式 SAVE 及冗余抵消 flags；记录三轮数值与耗时对照。 |
| 首版已实施并回归 | [rangular_openmp.md](implemented/rangular_openmp.md) | 独立 rangular_omp、线程私有角动量工作区、有界结果块与有序提交；1/2/4/8 线程及回退路径字节一致，23,459-CSF 单次 8 线程约 2.76×；排序/输出仍串行，非全套 GRASP 或物理验收。 |
| 已实施，运行器回归通过 | [rangular_version_comparison.md](implemented/rangular_version_comparison.md) | 指定原版/新版可执行文件，独立目录运行文本及可选描述符路径；保存日志、墙钟、阶段计时与 MCP 对照；7 项维护测试通过，支持差异/超时/提前退出检测。 |
| 已实施，输入/输出回归通过 | [rangular_descriptor_input.md](implemented/rangular_descriptor_input.md) | suite 预处理二进制、Fortran reader、文本数组对照和串行阶段计时；14 项专项测试、原生 rCSFs 171 CSF 对照及小型 MCP 一致性通过；生产性能、完整科学结果与 OpenMP 未验证。 |
| 首版输入和 OpenMP 已实施 | [rangular_descriptor_openmp_2026-11-10.md](plans/rangular_descriptor_openmp_2026-11-10.md) | 保留设计历史；输入、线程私有工作区、动态任务与有序提交已有首版，显式 context 和排序并行仍待后续工作。 |
| 已实施，解析器回归通过 | [mcp_compare.md](implemented/mcp_compare.md) | 2026-11-10，普通 rangular 串行/MPI MCP 解析、结构校验、按 CSF 对比较和 JSON 导出；合成回归与现有文件读回通过，不代表候选 rangular 或物理计算验收。 |
| 现行流程，静态核查完成 | [rangular_flow_and_state.md](implemented/rangular_flow_and_state.md) | 2026-11-10，基准 d181bf4；普通串行 rangular 的 block/CSF 对遍历、角系数缓冲与排序、共享工作区和文件状态。未实施 OpenMP，未运行性能或物理验证。 |
| 现行流程，静态核查完成 | [rmcdhf_flow_and_state.md](implemented/rmcdhf_flow_and_state.md) | 2026-11-10，基准 d181bf4；普通串行 rmcdhf 的矩阵/混合系数、径向积分、轨道改善、停止条件及共享状态；区分 strict 行为并标出 EAL 定义路径待复核。未实施 OpenMP，未运行性能或物理验证。 |
| 初步分析完成，诊断待实施 | [rmcdhf_mpi_single_node_diagnostics_2026-11-09.md](plans/rmcdhf_mpi_single_node_diagnostics_2026-11-09.md) | 普通/memory MPI 在大规模 CSF 下全程反复 D 状态；记录单节点条件、SETLAG 冻结侧冗余和 Davidson 候选。根因未确定，自动采样、服务器对照及 MPI/OpenMP 优化均未实施。 |
| 优先警告修复第一批通过 | [f2018_warning_repairs_2026-11-09.md](reviews/f2018_warning_repairs_2026-11-09.md) | GETALDWT、RCI 计数、Davidson、DACON、INTERP 及不可达步长分支；正确性复核 515→475，INTENT 9→0；32 配置、1,960 次 CTest 全通过，全开 93/93。 |
| 首次矩阵/基线设施证据 | [f2018_matrix_warning_baseline_2026-11-09.md](reviews/f2018_matrix_warning_baseline_2026-11-09.md) | 保留初次 32 配置、1,560 次 CTest 通过与 3,040 文件/515 条正确性诊断记录；当前基线及修复见上一项。 |
| 已更新，剩余正确性复核待完成 | [gfortran15-warning-baseline.md](reference/f2018/gfortran15-warning-baseline.md) | 当前 3,046 文件、4,804 条诊断，475 条正确性复核项；含无警告文件、配置次数和 source hash，不把登记视为已解决。 |
| strict 修复阶段验收通过 | [strict_scnsty_normalization_2026-11-09.md](reviews/strict_scnsty_normalization_2026-11-09.md) | SCNSTY 改为归一化后、阻尼前比较；公式/阈值/两轮条件不变，1/2-rank 均 11 轮收敛；保留当时 74/74 和旧新数值证据。 |
| Linux 初始修复证据 | [f2018_linux_davidson_restart_2026-11-09.md](reviews/f2018_linux_davidson_restart_2026-11-09.md) | 保留 Davidson、空 rank、memory IOS 修复及 71/73、strict 旧平台重现证据；后续 strict 时序修复及验收见上一项。 |
| 此前严格接口验证 | [f2018_strict_interfaces_npol_2026-11-09.md](reviews/f2018_strict_interfaces_npol_2026-11-09.md) | 内部接口/NP 路径/SOLVE 意图修复；保留 macOS Check/Release 各 64/68、4 skip 及一轮全输出一致记录；后续 Linux 状态见上一项。 |
| 设施初始验证记录 | [f2018_acceptance_facilities_2026-11-09.md](reviews/f2018_acceptance_facilities_2026-11-09.md) | profiles、强制清单、模块/安装隔离与 DDOT 14/14；保留初始 audit 50/63 及严格编译失败证据，后续状态见上一项。 |
| 已移除参数放行 | [argument_interfaces_2026-11-08.md](reviews/argument_interfaces_2026-11-08.md) | gfortran 15/F2018 全量构建不再使用 fallow-argument-mismatch；MPI/BLAS 显式接口、类型和有界数组区段修复；当前 CTest 36/47，11 项为既有失败；orbopt 已退役。 |
| 已移除旧程序 | [orbopt_retirement.md](implemented/orbopt_retirement.md) | 2026-11-08 删除四个 rmcdhf_orbopt 目标/源码和专用驱动，保留原版、strict/diag、CSFG；共有测试迁到 test/rmcdhf，suite 网格入口调整为 21 个。 |
| 已实施兼容性子批次 | [f2018_compatibility.md](implemented/f2018_compatibility.md) | 首轮标准语法及后续显式接口、NP 路径修复；严格 Check/Release 全量编译/链接通过，完整数值及静态状态验收未完成。 |
| 已完成本轮验证 | [f2018_gfortran15_batches_2026-11-08.md](reviews/f2018_gfortran15_batches_2026-11-08.md) | 全量构建通过；五项新增回归通过；四个 solver 一轮 SCF 共 28 文件与旧版本一致；CTest 33/44，11 个失败与旧基线相同，无新增失败。 |
| 已完成首轮编译诊断 | [f2018_gfortran15_compile_probe_2026-11-08.md](reviews/f2018_gfortran15_compile_probe_2026-11-08.md) | GNU 构建已加 `-std=f2018`，全新全量构建失败；15 个文件直接报错。补充生产源码诊断为 59 个文件失败，移除参数放行后为 143 个；列出共享库、平台函数、kind/语法与 MPI/BLAS 接口优先项，未修改计算实现或运行数值测试。 |
| 已实施 | [standard_grasp_integration.md](implemented/standard_grasp_integration.md) | 2026-11-08 导入普通 rangular、rwfnestimate、rbiotransform、rci、rtransition 五类程序共十个目标，复用 suite 共享库/参数；移除初猜程序局部容量定义，补丁器覆盖 25 个网格/4 个续算入口。完整 F2018 与原子物理验收尚未完成。 |
| 已完成接入检查 | [standard_grasp_integration_validation_2026-11-08.md](reviews/standard_grasp_integration_validation_2026-11-08.md) | GNU 15.3.0 全量及 serial-only 构建通过；35 项补丁器测试通过；CTest 28/39 通过，11 项 NP/strict 失败在未修改基线复现。记录既有并行模块输出冲突；不宣称全绿或物理验收。 |
| 已实施 | [jj2lsj_postprocessing.md](implemented/jj2lsj_postprocessing.md) | 2026-11-04 分阶段计时、稳定归并排序及 ASF 排序并行；唯一标签和输出按序执行，跨 block 流水线经计时评估暂缓。 |
| 已实施 | [jj2lsj_openmp.md](implemented/jj2lsj_openmp.md) | 2026-11-02 `jj2lsj_omp` 普通 CSF 的 OpenMP LS 行转换，保留默认提前停止与串行 `.lsj.T`；数值/文件回归及 Debug、无 OpenMP 构建通过；Gd I 真实输入兼容性通过，已有单次性能观测。 |
| 已实施 | [graspg_integration.md](implemented/graspg_integration.md) | 2026-10-30 GRASPG 计算源码、16 个程序、共享参数与串行/MPI CMake；工程回归通过，真实原子 SCF/RCI 验证未执行。 |
| 已实施，小型严格重启通过 | [rmcdhf_mem_strict_mpi.md](implemented/rmcdhf_mem_strict_mpi.md) | strict 比较时序修复后，O0/O2、1/2-rank 均在 11 轮满足原两轮判据；单轮/重启读写及解析反例通过。生产物理验证未完成。 |
| 已实施 | [rmcdhf_mem_diag_mpi.md](implemented/rmcdhf_mem_diag_mpi.md) | 2026-10-13 从原版独立复制的 P0/START 诊断目标；记录起始解、内部迭代和耦合系数，保留原算法及停止判据。START 逐值对照、投影边界及 Linux 1/2-rank 输出不变性回归通过；Ni 失稳根因尚待诊断重放。 |
| 已实施 | [rdensity_mpi.md](implemented/rdensity_mpi.md) | 2026-10-09 新增串行/MPI 自然轨道程序；CSF 分工、矩阵归约、解析小模型及 Ni 对照测试。大规模性能与物理收敛仍待验证。 |
| 已实施 | [common_parameters.md](implemented/common_parameters.md) | 2026-10-13 保留参数转导层，恢复原版参数值；2026-11-08 退役 orbopt 后现有 21 个网格入口、4 个 RCI 续算入口，NNNP=590、NNN1=600。 |
| 已实施 | [RHFS_gJ_report.md](implemented/RHFS_gJ_report.md) | 串行 RHFS 的加载、矩阵元、ASF 投影及输出链仍存在。仅覆盖该路径，不是 gj90/MPI 完整说明。 |
| 已实施 | [rmcdhf_nuclear_polarization.md](implemented/rmcdhf_nuclear_polarization.md) | 2026-10-24 从 grasp 新复制四个原版 rmcdhf 目标并接入局域 E1/E2 NP；已有 orbopt 保持原状：独立输入、共享势、root I/O 与集体失败协议；串行/MPI 代码测试 9/9、14/14 通过，物理计算验证未执行。 |
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
- 二十五个入口（含普通 GRASP 计算链、四个 GRASPG、新复制的四个原版 RMCDHF 及 P0 诊断目标）保留各自的原版点核/有限核默认值与 ACCY 初始化时序；
  当前有限核默认 `N=NNNP=590`。suite 脚本布局按需修改集中容量配置和这些入口，保持转导接口。
- [test/CMakeLists.txt](../test/CMakeLists.txt) 保留共享容量检查和网格脚本回归；
  历史测试数量不能代替当前实际验证。严格内存 MPI 的 Linux 重启测试另见其说明。
- 历史 `N=2990,H=0.05` 的风险审查与候选网格实验仍有参考价值，
  但该组合已不是仓库默认值，参数恢复本身也不代表生产体系收敛。

### RHFS 与独立 gj90

- [hfs92.f90](../src/appl/rhfs90/hfs92.f90) 与
  [hfsgg.f90](../src/appl/rhfs90/hfsgg.f90) 对应现行串行流程报告。
- [gj90.f90](../src/appl/gj90/gj90.f90)、
  [compute_gj_only.f90](../src/appl/gj90/compute_gj_only.f90) 和
  [setout_gj.f90](../src/appl/gj90/setout_gj.f90) 已实现独立入口、NVEC 对角累加
  和 `.gj/.cgj` 输出，取代蓝图的“待新建”描述。
- [compute_gj_only_mpi.f90](../src/appl/gj90/compute_gj_only_mpi.f90) 已按
  `IC=MYID+1,NCF,NPROCS` 分工并归约；主循环 OpenMP/线程工作区重构未见实现。
- 蓝图曾把 RHFS `getmixblock.f90` 写在 `lib9290/` 下，实际文件位于
  [src/appl/rhfs90](../src/appl/rhfs90/getmixblock.f90)，归档正文已更正路径。

### RHFS/RIS 性能优化

- [hfsgg_mpi.f90](../src/appl/rhfs90/hfsgg_mpi.f90) 保留三角/宇称、内联 IQA
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
2. 当前 [CMakeLists.txt](../CMakeLists.txt) 的 gj90 数值回归注册依赖
   `tests/gj90_regression.sh` 和 `tests/gj90_mpi_regression.sh`，但当前受版本控制
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
