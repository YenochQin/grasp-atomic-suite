# 文档索引与实施状态审计

核查日期：2026-10-05。源码及原文基准：`8b298a6`。
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
├── plans/          # 仍需执行的实验/验证方案
├── reviews/        # 已完成的审查与证据记录
├── reference/      # 理论背景和参考草稿
└── archive/        # 已被替代的蓝图、迁移快照和阶段性报告
```

当前运行与构建先看[仓库 README](../README_ZH.md)和
[脚本说明](../scripts/README.md)；参数管理看 `implemented/`，
实施网格实验看 `plans/`，追溯变更和旧测量看 `archive/`。
归档文档保留设计与实验历史，不能直接当成当前配置指南。

## 逐份文档结论

| 分类 | 文档 | 实施状态与适用范围 |
| --- | --- | --- |
| 已实施 | [common_parameters.md](implemented/common_parameters.md) | 参数集中管理和六个入口已接入。当前默认值有效；历史构建/测试记录不等于物理收敛。 |
| 已实施 | [RHFS_gJ_report.md](implemented/RHFS_gJ_report.md) | 串行 RHFS 的加载、矩阵元、ASF 投影及输出链仍存在。仅覆盖该路径，不是 gj90/MPI 完整说明。 |
| 待执行 | [grasp_grid_parameter_tuning.md](plans/grasp_grid_parameter_tuning.md) | 修改脚本已实现，文档配置已核查；D/L/O/ACCY 和生产模型物理验收仍待执行。 |
| 已完成审查 | [grasp_grid_methodology_review.md](reviews/grasp_grid_methodology_review.md) | 源码审查和最小例程诊断已完成，文档建议已落实。旧 P 配置等问题描述有指定历史版本。 |
| 理论参考 | [lande_g.md](reference/lande_g.md) | 概念、ASF 与 LS 近似的参考草稿，不能划为功能实施计划。算符规范与系数尚未完整独立核对。 |
| 历史蓝图 | [gj90_blueprint.md](archive/gj90_blueprint.md) | 阶段 A/C 主体已实现，B 有历史样例对比记录；阶段 D 的 OpenMP/线程重构未见实现。不再是当前开发清单。 |
| 历史阶段报告 | [RHFS_MPI_hfsgg_optimization.md](archive/RHFS_MPI_hfsgg_optimization.md) | 前置筛选仍在；诊断计数、直接 ICHKQ1 调用和“当前瓶颈”描述已被后续代码覆盖。 |
| 历史分支总结 | [1.0.1-dev.1_optimization_summary.md](archive/1.0.1-dev.1_optimization_summary.md) | HFSGG 与 DIAGA3/4 保留。RECOP2 独立快速预检查未保留，已修正原文误记；旧计时未复测。 |
| 历史分支总结 | [1.1.1-dev.1_optimization_summary.md](archive/1.1.1-dev.1_optimization_summary.md) | RIS 主要优化仍存在；历史耗时、测试名称和测试覆盖不代表当前数值验收。 |
| 历史迁移记录 | [rmcdhf_migration.md](archive/rmcdhf_migration.md) | 四个 orbopt 目标和共享库已落地；性质程序 590 点默认、多入口补丁及旧测试数量已被集中管理替代。 |

## 核查依据与具体过时项

### 参数与仓库迁移

- [suite_parameters_M.f90](../src/lib/libmod/suite_parameters_M.f90) 是公共参数定义；
  [parameter_def_M.f90](../src/lib/libmod/parameter_def_M.f90) 只转导容量名。
- [radial_grid_defaults_M.f90](../src/lib/libmod/radial_grid_defaults_M.f90) 被四个
  RMCDHF、RHFS 和 RIS 网格入口调用；g_J 复用 RHFS 的入口。
- 当前全部有限核默认 `N=NNNP=2990`。迁移记录里的性质程序
  `N=MIN(590,NNNP)` 仅适用于早期状态；当前脚本的 suite 布局只修改集中参数文件。
- [test/CMakeLists.txt](../test/CMakeLists.txt) 已加入共享容量与默认值逻辑检查。
  迁移记录的 6/5 项、12 个离线用例和公共参数说明的 16 个用例都是历史数量，
  当前登记与实际覆盖须查看 CMake 和测试源码，不能从旧报告推断。
- 默认 `N=2990,H=0.05` 的数值风险与候选组见网格实验方案；代码已经统一，
  不代表这一默认组合已通过物理验证。

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
[stage guard 实施记录](../test/rmcdhf_orbopt/STAGE_GUARD_IMPLEMENTATION.md)，
没有把大型外部计算记录迁入 `docs/`。

## 后续维护规则

新增文档先按用途选择这五类目录，写明日期、代码版本、已实施部分和待验证部分。
完成计划时更新状态并指向实现证据；归档时保留历史参数和测量背景，注明当前替代入口。
源码变更触及现行说明时同步更新文档与此索引。仓库内部链接使用相对路径，
避免依赖某台机器的绝对目录和已经漂移的源码行号。
