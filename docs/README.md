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

## 逐份文档结论

| 分类 | 文档 | 实施状态与适用范围 |
| --- | --- | --- |
| 已实施 | [rmcdhf_mem_strict_mpi.md](implemented/rmcdhf_mem_strict_mpi.md) | 2026-10-13 新增独立内存 MPI 严格停止目标；编译、判据单元测试通过。完整 Linux MPI 重启与生产物理验证待执行。 |
| 已实施 | [rmcdhf_mem_diag_mpi.md](implemented/rmcdhf_mem_diag_mpi.md) | 2026-10-13 从原版独立复制的 P0/START 诊断目标；记录起始解、内部迭代和耦合系数，保留原算法及停止判据。START 逐值对照、投影边界及 Linux 1/2-rank 输出不变性回归通过；Ni 失稳根因尚待诊断重放。 |
| 已实施 | [rdensity_mpi.md](implemented/rdensity_mpi.md) | 2026-10-09 新增串行/MPI 自然轨道程序；CSF 分工、矩阵归约、解析小模型及 Ni 对照测试。大规模性能与物理收敛仍待验证。 |
| 已实施 | [common_parameters.md](implemented/common_parameters.md) | 2026-10-13 保留参数转导层，恢复原版参数值；2026-10-24 新增四个原版 RMCDHF 入口后共十三个入口，NNNP=590、NNN1=600。 |
| 已实施 | [RHFS_gJ_report.md](implemented/RHFS_gJ_report.md) | 串行 RHFS 的加载、矩阵元、ASF 投影及输出链仍存在。仅覆盖该路径，不是 gj90/MPI 完整说明。 |
| 已实施 | [rmcdhf_nuclear_polarization.md](implemented/rmcdhf_nuclear_polarization.md) | 2026-10-24 从 grasp 新复制四个原版 rmcdhf 目标并接入局域 E1/E2 NP；已有 orbopt 保持原状：独立输入、共享势、root I/O 与集体失败协议；串行/MPI 代码测试 9/9、14/14 通过，物理计算验证未执行。 |
| 已完成代码计划 | [rmcdhf_nuclear_polarization_implementation.md](plans/rmcdhf_nuclear_polarization_implementation.md) | 保留实施设计、完成清单与结果；现行输入契约见已实施说明。 |
| 待执行 | [grasp_grid_parameter_tuning.md](plans/grasp_grid_parameter_tuning.md) | 修改脚本已实现，文档配置已核查；D/L/O/ACCY 和生产模型物理验收仍待执行。 |
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
- 十三个入口（含新复制的四个原版 RMCDHF 及 P0 诊断目标）保留各自的原版点核/有限核默认值与 ACCY 初始化时序；
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

`test/rmcdhf_orbopt/` 中的测试说明和实验结果保留在测试旁，
分类索引仅链接[测试入口](../test/rmcdhf_orbopt/README.md)、
[历史结果](../test/rmcdhf_orbopt/RESULTS.md)及
[已退役的 stage guard 历史记录](archive/orbopt_stage_guard_retired.md)，
2026-10-19 移除 Python 阶段验收流程及其专属入口，实施记录已归档；
Fortran 求解器、构建目标和直接回归流程保留。没有迁入大型外部计算数据。

## 后续维护规则

新增文档先按用途选择这五类目录，写明日期、代码版本、已实施部分和待验证部分。
完成计划时更新状态并指向实现证据；归档时保留历史参数和测量背景，注明当前替代入口。
源码变更触及现行说明时同步更新文档与此索引。仓库内部链接使用相对路径，
避免依赖某台机器的绝对目录和已经漂移的源码行号。
