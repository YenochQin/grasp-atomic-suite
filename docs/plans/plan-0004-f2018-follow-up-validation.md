# F2018 后续数值与状态验收

以 Suite `0f35673` 的保留目标为范围。接口模块化、`USE ONLY`、数字标签 FORMAT 内联、
首批变量语义迁移和 V3 输入替换均有实施记录。构建和告警验收使用 GNU Fortran 16.2。
[构建策略](../implemented/design-0002-fortran-build-profiles.md) · [已实施兼容性](../implemented/design-0004-f2018-compatibility.md)

| 剩余问题 | 验证范围 |
| --- | --- |
| Fortran 接口、状态和诊断 | 以现行源码重扫 fortls/Fortitude 和编译诊断。逐项核查 `INTENT`、可见性、公开可变状态、回调、数组边界、实际重导出、固定单元号和错误传播。按调用链核实跨调用局部状态、MPI 计数与绑定边界。 |
| 串行与 MPI 的多根数值 | 用同构建参数、同基组及初始轨道，比较各保留程序与自身版本基线。记录多根权重与态匹配、每轮混合系数及能量、矩阵或本征残差、轨道正交性和收敛状态。涉及 MPI 时记录 rank 布局及归约容差。 |
| 完整性质链 | 固定轨道验证 RCI 的 Coulomb、Breit 与实际支持的质量/QED 分支。对需要覆盖的 RHFS、g_J、RIS、RDENSITY、双正交和跃迁路径分别核对目标性质及分支输入。RIS 缓存与单位的小回归已存在，完整性质对照仍待执行。 |
| 平台和数值假设 | 在实际部署平台核查编译器、MPI、BLAS/LAPACK 和运行结果。测试告警、运行时检查、数值回归及模型收敛各自记录。 |

Fortran 命名或诊断只按当前实体与源码核对。旧清单的 19,579 项和扫描统计对应当时源码。
GNU 16.2.1/Linux x86-64 复采的 37 个新增诊断身份中有 9 个 `-Wmaybe-uninitialized`，
涉及双正交轨道变换、RMCDHF 乘子、MCP 收缩数组边界和求解器索引。这里同时包含
未改源码上的平台差异和改动后位置变化。对每条诊断仍需核对定义路径与可达输入。
[当前告警索引](../reference/f2018/evidence-0019-gnu16-warning-baseline.md)

已核实的第一批 529 项可追溯到 [命名映射](../reference/non_csfg_variable_renames_2026-11-16.json)。
原扫描见 [Fortitude/fortls 报告](../reviews/evidence-0012-fortls-fortitude-audit.md)，
用例和容差协议见 [历史 GNU 15 迁移方案](../archive/plan-0002-f2018-migration-gnu15.md)。
