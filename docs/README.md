# 文档索引

2026-11-23 核对 Suite `0f35673` 的实现、验证记录和剩余计划。
“已实施”表示代码/接口存在，不表示所有体系的精度、收敛和性能已验收。

## 阅读入口

- 构建与运行：[仓库说明](../README_ZH.md)、[构建策略](implemented/design-0002-fortran-build-profiles.md)、
  [Rust 原生 CSF 输入](implemented/design-0001-native-csf-descriptors.md)、[脚本说明](../scripts/README.md)。
- 容量与网格：[共同参数](implemented/design-0003-shared-grid-parameters.md)；参数容量与运行时网格点数分别核定。
- 普通计算链：[rangular 流程](implemented/design-0017-rangular-flow-and-state.md) →
  [rmcdhf 流程](implemented/design-0019-rmcdhf-flow-and-state.md)。
- 后续验收：[径向网格](plans/plan-0001-radial-grid-validation.md)、[F2018 现代化实施计划](plans/plan-0004-f2018-follow-up-validation.md)、
  [MPI MCP 生产对照](plans/plan-0005-rmcdhf-mcp-production-io.md)。

`implemented/` 保留功能契约，`plans/` 保留未完成工作，`reviews/` 保留独立分析和验证证据，
`reference/` 保存理论及机器可读来源，`archive/` 保存仍有参考价值的历史背景。
旧输入说明和退役程序记录按其标注版本阅读，不能据此恢复旧生产入口。

## 现行功能说明与保留的实现背景

| 文档 | 范围与状态 |
| --- | --- |
| [design-0016-inline-output-formats.md](implemented/design-0016-inline-output-formats.md) | 数字标签 `FORMAT` 内联，覆盖共享库和应用的 903 处 `WRITE`、2 处 `PRINT`；源码、预处理和字节对照已验证。 |
| [design-0006-shared-library-naming.md](implemented/design-0006-shared-library-naming.md) | 共享状态与基础接口库的命名约定、Fortran 调用方迁移和完整模块对照。 |
| [design-0007-shared-library-procedure-names.md](implemented/design-0007-shared-library-procedure-names.md) | 共享库过程名称迁移表与现行 Fortran 调用入口。 |
| [design-0001-native-csf-descriptors.md](implemented/design-0001-native-csf-descriptors.md) | 现行 V3 生产 CSF 输入、ABI、MPI、双基组、低内存、错误传播和数值回归。接入设计与验收用例见 [已实施的 V3 设计与验收](plans/plan-0006-csf-descriptor-v3-integration.md)。 |
| [design-0018-rangular-angular-reuse.md](implemented/design-0018-rangular-angular-reuse.md) | 已实施，串行专项回归通过；默认关闭的有界 CSF 对缓存；保留顺序、占据/耦合及轨道共享关系，复用原始 T/V 系数，普通 MCP 写出不变；包含命中率、缓存容量与 RSS 测量，MPI 未接入；rangular OpenMP 已移除。 |
| [design-0023-rwfnestimate-header-input.md](implemented/design-0023-rwfnestimate-header-input.md) | header 输入说明；保留 header 轨道表、电子数和初猜检查背景；当前仅 header 输入，无旧文本回退。 |
| [adr-0002-retire-strict-and-diagnostic-rmcdhf.md](implemented/adr-0002-retire-strict-and-diagnostic-rmcdhf.md) | 已移除被否决变体；删除 strict/diag 两个独立求解器、专用测试、构建/安装入口及当前二进制；普通四种 RMCDHF 和共享库保持原流程。 |
| [design-0005-program-naming.md](implemented/design-0005-program-naming.md) | 普通程序命名；保留应用目录和可执行文件去掉 90 的命名契约；CSFG 已随后移除，共享库及上游来源命名保留。 |
| [design-0002-fortran-build-profiles.md](implemented/design-0002-fortran-build-profiles.md) | 现行构建策略；CSFG 已移除；生产默认 -O2，OpenMP 仅保留 jj2lsj；标准/接口/运行时诊断放 Check/Audit；serial、MPI 和 jj2lsj OpenMP 组成 24 个有效配置，16 个诊断配置沿用历史关闭 angular OpenMP 的预算 ID；历史预算范围不代表当前告警验收。 |
| [rangular 对照工具](../scripts/README.md) | 原生描述符默认和显式输入、历史文本 oracle、阶段计时与 MCP 对照。 |
| [design-0024-mcp-compare.md](implemented/design-0024-mcp-compare.md) | 已实施，解析器回归通过；2026-11-10，普通 rangular 串行/MPI MCP 解析、结构校验、按 CSF 对比较和 JSON 导出；合成回归与现有文件读回通过，不代表候选 rangular 或物理计算验收。 |
| [design-0017-rangular-flow-and-state.md](implemented/design-0017-rangular-flow-and-state.md) | 现行流程，静态核查完成；2026-11-10，基准 d181bf4；普通串行 rangular 的 block/CSF 对遍历、角系数缓冲与排序、共享工作区和文件状态。未实施 OpenMP，未运行性能或物理验证。 |
| [design-0019-rmcdhf-flow-and-state.md](implemented/design-0019-rmcdhf-flow-and-state.md) | 现行流程，静态核查完成；2026-11-10，基准 d181bf4；普通串行 rmcdhf 的矩阵/混合系数、径向积分、轨道改善、停止条件及共享状态；标出 EAL 定义路径待复核。未实施 OpenMP，未运行性能或物理验证。 |
| [design-0021-rmcdhf-mpi-diagnostics.md](implemented/design-0021-rmcdhf-mpi-diagnostics.md) | 外部采样首版已实现，服务器验证待进行；2026-11-15，普通 MPI 原版流程与兼容修复核查；Linux PID 树、D/R 线程、I/O、缺页、PSI 和 cgroup 只读采样。模拟回归通过；真实开销及停顿根因未验证。 |
| [adr-0001-retire-rmcdhf-orbopt.md](implemented/adr-0001-retire-rmcdhf-orbopt.md) | 已移除旧程序；2026-11-08 删除四个 rmcdhf_orbopt 目标/源码和专用驱动，记录 2026-11-08 的 orbopt 退役；strict/diag 后续也已移除，当前网格范围以共同参数说明为准。 |
| [design-0004-f2018-compatibility.md](implemented/design-0004-f2018-compatibility.md) | 现行兼容性与接口；汇总标准 kind、共享库模块 API、MPI 重导出例外、输入/数值修复、来源校验与验证限制。 |
| [design-0008-core-library-naming.md](implemented/design-0008-core-library-naming.md) | 公共基础库的目录、静态库、Fortran 文件与模块命名对照。 |
| [Core 计算公式](../src/lib/core/readme.md) | 公共基础库中的径向网格、求积、轨道、核势、角动量和矩阵运算公式。 |
| [design-0009-core-variable-migration.md](implemented/design-0009-core-variable-migration.md) | `libcore` 变量与常量的逐作用域对应表，包含形参、保留名称与共享状态别名。 |
| [design-0015-breit-angular-support-migration.md](implemented/design-0015-breit-angular-support-migration.md) | Breit 角系数支持库的命名迁移、共享径向积分系数过程和角系数收集接口。 |
| [design-0012-spin-angular-library-naming.md](implemented/design-0012-spin-angular-library-naming.md) | 自旋–角向库的目录、库目标、89 个源文件和 Fortran 模块命名对照。 |
| [design-0013-spin-angular-variable-migration.md](implemented/design-0013-spin-angular-variable-migration.md) | `spin_angular` 变量与常量的逐过程对照，涵盖子壳层位置、角动量编码、张量矩阵元和 Breit 角向系数。 |
| [design-0025-rlevels.md](implemented/design-0025-rlevels.md) | 原版 `rlevels` 的 Suite 接入，含 `.m`、`.cm`、LSJ 标签匹配规则和回归入口。 |
| [design-0029-jj2lsj-postprocessing.md](implemented/design-0029-jj2lsj-postprocessing.md) | 已实施；2026-11-04 分阶段计时、稳定归并排序及 ASF 排序并行；唯一标签和输出按序执行，跨 block 流水线经计时评估暂缓。 |
| [design-0028-jj2lsj-openmp.md](implemented/design-0028-jj2lsj-openmp.md) | 已实施；2026-11-02 `jj2lsj_omp` 普通 CSF 的 OpenMP LS 行转换，保留默认提前停止与串行 `.lsj.T`；数值/文件回归及 Debug、无 OpenMP 构建通过；Gd I 真实输入兼容性通过，已有单次性能观测。 |
| [design-0026-rdensity-mpi.md](implemented/design-0026-rdensity-mpi.md) | 已实施；2026-10-09 新增串行/MPI 自然轨道程序；CSF 分工、矩阵归约、解析小模型及 Ni 对照测试。大规模性能与物理收敛仍待验证。 |
| [design-0003-shared-grid-parameters.md](implemented/design-0003-shared-grid-parameters.md) | 已实施；2026-11-16 使用语义容量名称，由 cal_constants 经 grasp_parameters 转导；CSFG 移除后覆盖 15 个普通网格入口、2 个 RCI 续算入口，默认容量仍为 590/600。 |
| [design-0014-davidson-eigensolver-migration.md](implemented/design-0014-davidson-eigensolver-migration.md) | Davidson 求解库的目录、静态库、Fortran 文件、模块和变量命名对照。 |
| [design-0010-mpi-library-naming.md](implemented/design-0010-mpi-library-naming.md) | MPI 支持库的目录、静态库、Fortran 文件和模块命名对照。 |
| [design-0011-mpi-variable-migration.md](implemented/design-0011-mpi-variable-migration.md) | MPI 支持库的作用域级变量与常量对应表，含通信、输入、目录管理和稀疏矩阵例程。 |
| [design-0027-rhfs-lande-g.md](implemented/design-0027-rhfs-lande-g.md) | 已实施；串行 RHFS 的加载、矩阵元、ASF 投影及输出链仍存在。仅覆盖该路径，不是 gj/MPI 完整说明。 |
| [design-0022-rmcdhf-nuclear-polarization.md](implemented/design-0022-rmcdhf-nuclear-polarization.md) | 已实施，普通 MPI 已移除 NP 入口；2026-11-15 普通 rmcdhf_mpi 恢复直接核势初始化；rmcdhf、rmcdhf_mem、rmcdhf_mem_mpi 保留可选 NP。保留 2026-10-24 历史代码验证，物理计算验证未执行。 |
| [design-0020-rmcdhf-mcp-io.md](implemented/design-0020-rmcdhf-mcp-io.md) | MPI 计算链统一使用目录格式。普通版使用每轮 EOL 收缩，内存版读取相同分片。生产 A/C 的 SCF 数值和墙钟对照已完成。 |

## 未完成的计划与清单

| 文档 | 范围与状态 |
| --- | --- |
| [RMCDHF 的 MCP 读取优化实施计划](plans/plan-0005-rmcdhf-mcp-production-io.md) | 阶段 0–2 已通过本地回归。阶段 3 已完成 Zn I 服务器矩阵分段和逐轮驻留诊断、Cd I 首组冷态 A/C 对照。Cd I 重复配对、生成期磁盘峰值、资源预算及 B 对照待验收。 |
| [plan-0003-rmcdhf-mpi-stall-diagnostics.md](plans/plan-0003-rmcdhf-mpi-stall-diagnostics.md) | 诊断计划，D0 首版已实现；普通/memory MPI 在大规模 CSF 下全程反复 D 状态；记录单节点条件、SETLAG 冻结侧冗余和 Davidson 候选。普通 MPI 已有 MCP 阶段计时。服务器联合采样、扩展性与根因诊断待完成。 |
| [plan-0001-radial-grid-validation.md](plans/plan-0001-radial-grid-validation.md) | 待执行；修改脚本已实现，文档配置已核查；D/L/O/ACCY 和生产模型物理验收仍待执行。 |
| [F2018 现代化实施与验收计划](plans/plan-0004-f2018-follow-up-validation.md) | P0 续算与 V3 基线、P1 MPI 绑定及 BLAS/LAPACK LP64 配置验收已完成。四种 RMCDHF 已取得 Zn I AS1 的七根收敛轨迹及逐 rank 峰值内存，固定轨道 RCI 已核对 Ni I 双根的稀疏内存与磁盘结果。双正交已核对单电子轨道、态向量和角缓存；性质与跃迁及候选版本完整矩阵待验收。 |

## 实施过程记录

- [2026-11-23 RMCDHF MCP 阶段 3 实施过程](plans/2026-11-23-plan-0005-rmcdhf-mcp-phase3.md)记录基准来源、采样包演进、复跑入口和证据位置。

## 分析与验证证据

| 文档 | 范围与状态 |
| --- | --- |
| [双正交轨道、态向量和角缓存](reviews/evidence-0027-biorthogonal-states-and-cache.md) | 串行与 2-rank MPI 的正逆序及不等维解析对照，错误基组缓存的波函数偏差和重算验证。 |
| [Zn I 多根 RMCDHF 收敛与逐 rank 内存](reviews/evidence-0026-rmcdhf-zn-multiroot.md) | 四种程序的七根能量、混合系数、径向轨迹及串行／2-rank 峰值 RSS 对照。 |
| [RMCDHF MCP 阶段 3 实测与续跑](reviews/evidence-0025-rmcdhf-mcp-phase3-benchmark.md) | Zn I AS2 在服务器 XFS 和本机同盘 Btrfs、XFS 的 A/C 实测，以及页驻留诊断与外部 Rust 静态库复用。 |
| [多 J 块 MCP 容量与读写成本](reviews/evidence-0017-rangular-mcp-io.md) | 串行与双进程文件字节对照，48-rank 历史作业的容量和读写核算，以及面向 RMCDHF 的分块积分目录、每次密度更新后的流式收缩和验收方法。 |
| [evidence-0016-fortran-source-and-readme-review.md](reviews/evidence-0016-fortran-source-and-readme-review.md) | 已修复审查列出的 Fortran 命名、注释及 README 问题，记录串行和 MPI 回归及 macOS 内存探测限制。 |
| [evidence-0015-fortran-alias-audit.md](reviews/evidence-0015-fortran-alias-audit.md) | Fortran 导入重命名与旧命名兼容入口已清除，含迁移对照、旧接口拒绝验证及运行结果比较。 |
| [evidence-0014-native-csf-descriptor-validation.md](reviews/evidence-0014-native-csf-descriptor-validation.md) | 原生 CSF 输入的功能回归、矩阵及告警门禁边界；保留独立验收证据。 |
| [evidence-0012-fortls-fortitude-audit.md](reviews/evidence-0012-fortls-fortitude-audit.md) | 迁移前工具审查；基准 5806884；保留原静态扫描及解析边界。之后生产 CSF reader 已替换，不能将旧风险与计数当作当前状态。 |
| [evidence-0013-fortls-fortitude-rescan.md](reviews/evidence-0013-fortls-fortitude-rescan.md) | 工具对照证据；保留 S1–S4 后的解析器/编译器对照与检查边界；旧计数和 CSF 输入风险已由后续修改覆盖，不能直接作为当前执行清单。 |
| [evidence-0011-rmcdhf-mcp-precontraction.md](reviews/evidence-0011-rmcdhf-mcp-precontraction.md) | 已完成依赖与预收缩分析；普通 RMCDHF 的矩阵/势两条 MCP 消费路径；EOL 动态密度限制、全冻结积分收缩为 H_fixed、角向模板和迭代期收缩缓存。数学可行性与现行支持分开，未改代码或运行收缩/SCF 实验。 |
| [evidence-0010-rangular-mpi-mcp-storage.md](reviews/evidence-0010-rangular-mpi-mcp-storage.md) | 已完成源码分析与文件压缩实验；普通 MCP 的 rank/通道拆分、记录布局、临时峰值及 RMCDHF 读取约束；已有串行样本 Zstd 无损缩减 49.6%–64.5%，逐文件 SHA-256 一致。容器、分块压缩和冗余列号消除均未实施，生产 MPI/SCF 性能待测。 |
| [evidence-0009-rangular-reuse-benchmark.md](reviews/evidence-0009-rangular-reuse-benchmark.md) | 已完成输出与性能对照；独立构建修改前基准，两组数据各三轮，文本/描述符 MCP 零容差及字节一致；合成样本命中 99.93%、整程 1.21×，现有 14,551-CSF 样本几乎不命中、慢约 37%；默认保持关闭。 |
| [evidence-0008-rwfnestimate-nan-diagnosis.md](reviews/evidence-0008-rwfnestimate-nan-diagnosis.md) | 根因已复现并修复，物理优化待验证；Ni/Zn 首次 h/i 的搜索上界变负，导致 SETPOT/TAIL 警告及 NaN 误报成功；记录数值链、suite 迁移、上游恢复；同步命名迁移及退役更新后，Release 全量 66/66 通过。 |
| [evidence-0005-rangular-mpi-openmp-benchmark.md](reviews/evidence-0005-rangular-mpi-openmp-benchmark.md) | 已完成单机性能与输出对照；基于 c5f36d2，M4 上 14,551 CSFs、33 次 Release 运行；8 并发 MPI 4.06 秒、OpenMP 10.56 秒；描述符主要加快输入，整程收益较小；串行/线程字节一致，MPI 零容差等价。 |
| [evidence-0006-rangular-mpi-openmp-x86-diagnosis.md](reviews/evidence-0006-rangular-mpi-openmp-x86-diagnosis.md) | 已完成 x86 性能诊断；基于 c47f4d4，23,459 CSFs 三轮：OMP 8 线程中位数 27.06 秒、MPI 8 进程 13.49 秒；38,461,197 项贡献零差异。隔离阶段探针确认并行范围、逐列批次等待和串行结果处理的限制；未修改生产计算代码。 |
| [evidence-0004-rangular-mpi-descriptor-benchmark.md](reviews/evidence-0004-rangular-mpi-descriptor-benchmark.md) | 已完成 MPI 描述符测量；同一 8-rank 程序、23,459 CSFs 复核三轮：输入含广播中位数 0.411→0.061 秒，整程 12.565→12.162 秒，约减少 3.2%；96 个 MCP 文件与原 MPI 输出一致。 |
| [evidence-0007-rmcdhf-memory-variant-comparison.md](reviews/evidence-0007-rmcdhf-memory-variant-comparison.md) | 已完成源码对比，数值影响待验证；基准 db39505；四种普通 RMCDHF 的 MCP 缓存、乘子轨道对筛选、MPI 乘子公式和径向边界差异。不能视为纯存储切换；未运行原子计算或性能基准。 |
| [evidence-0003-jj2lsj-postprocessing-validation.md](reviews/evidence-0003-jj2lsj-postprocessing-validation.md) | 已完成验证；AS4 稳定排序及整理约 9.60→0.045 s，总耗时单次约减少 2.1%；AS4/AS1 输出一致，标签异步化经计时评估暂缓。 |
| [evidence-0002-jj2lsj-gd-i-validation.md](reviews/evidence-0002-jj2lsj-gd-i-validation.md) | 已完成验证；六套 Gd I 输入的 18 个输出文件一致；两个代表案例的 994,894 个 LS 系数在 4/8 线程下与串行逐位一致，包含真实提前停止路径。 |
| [evidence-0001-radial-grid-methodology-review.md](reviews/evidence-0001-radial-grid-methodology-review.md) | 已完成审查；源码审查和最小例程诊断已完成，文档建议已落实。旧 P 配置等问题描述有指定历史版本。 |

## 参考资料

| 文档 | 范围与状态 |
| --- | --- |
| [evidence-0019-gnu16-warning-baseline.md](reference/f2018/evidence-0019-gnu16-warning-baseline.md) | GNU Fortran 16.2.1 Linux x86-64 复采，覆盖 16 个普通 Check/Audit 配置及带目录 MCP 目标。逐条正确性复核仍在 F2018 后续验收范围内。 |
| [evidence-0018-gnu15-warning-baseline.md](reference/f2018/evidence-0018-gnu15-warning-baseline.md) | GNU Fortran 15.3.0 的历史告警记录，包含当时 24 个配置的编译覆盖与诊断。 |
| [design-0032-lande-g-theory.md](reference/design-0032-lande-g-theory.md) | 理论参考；概念、ASF 与 LS 近似的参考草稿，不能划为功能实施计划。算符规范与系数尚未完整独立核对。 |
| [design-0033-nuclear-polarization-parameters.md](reference/design-0033-nuclear-polarization-parameters.md) | 参数参考；2026-10-24 核对 Flambaum 2021 的 E1/E2 强度、50 keV 修正极化率、截断参数与原子单位换算；记录拟合系数内部差异，未建立生产核素参数库或物理精度结论。 |

## 历史背景

| 文档 | 范围与状态 |
| --- | --- |
| [普通 GRASP 导入](archive/2026-11-08-plan-0002-standard-grasp-integration.md) | 保存原始程序来源、接入范围及历史构建记录。 |
| [GNU 15 迁移方案](archive/plan-0002-f2018-migration-gnu15.md) | 保留数值夹具与验收协议，旧构建和目标按原版本阅读。 |
| [jj2lsj 后处理方案](archive/design-0030-jj2lsj-postprocessing-pipeline.md) | 稳定排序、计时和 ASF 排序已实施，跨块流水线经计时后暂缓。 |
| [design-0031-gj-program-blueprint.md](archive/design-0031-gj-program-blueprint.md) | 历史蓝图；阶段 A/C 主体已实现，B 有历史样例对比记录；阶段 D 的 OpenMP/线程重构未见实现。不再是当前开发清单。 |
| [evidence-0023-rhfs-mpi-hfsgg-optimization.md](archive/evidence-0023-rhfs-mpi-hfsgg-optimization.md) | 历史阶段报告；前置筛选仍在；诊断计数、直接 ICHKQ1 调用和“当前瓶颈”描述已被后续代码覆盖。 |
| [evidence-0021-v1-0-1-dev-1-optimization-summary.md](archive/evidence-0021-v1-0-1-dev-1-optimization-summary.md) | 历史分支总结；HFSGG 与 DIAGA3/4 保留。RECOP2 独立快速预检查未保留，已修正原文误记；旧计时未复测。 |
| [evidence-0022-v1-1-1-dev-1-optimization-summary.md](archive/evidence-0022-v1-1-1-dev-1-optimization-summary.md) | 历史分支总结；RIS 主要优化仍存在；历史耗时、测试名称和测试覆盖不代表当前数值验收。 |

补充历史资料包括 [orbopt 实验结果](archive/orbopt/evidence-0024-orbopt-as2-diagnostics.md)、
[GRASPG 上游说明](reference/graspg/evidence-0020-graspg-upstream-readme.md)
和 [GRASPG 源码清单](reference/graspg/source-manifest.json)。

## 验证边界与维护

网格及活性空间收敛、真实多根/完整性质计算、RIS 数值对照、服务器 MPI 性能及理论算符约定，
按各自计划与说明中的未完成范围继续核定。历史测试数不代表当前自动回归覆盖。
现有 gj 回归脚本缺口及 RIS 专项数值覆盖问题不能靠旧优化报告关闭；测试入口见
[test/CMakeLists.txt](../test/CMakeLists.txt)和[求解器测试说明](../test/rmcdhf/README.md)。

已完成的 USE ONLY、格式整理、大小写/首批改名、模块化、输入/数值修复和 F2018 分批报告，
以及完成的 CSF/NP 代码计划不再逐批保留；实现约束集中在
[兼容性说明](implemented/design-0004-f2018-compatibility.md)、原生输入及 NP 说明，原始记录由 Git 历史追溯。
测试引用的来源快照、告警基线和仍供后续工作使用的清单继续保留。

后续计划完成时，将仍适用的接口、操作方式及验证限制并入现行说明，删除重复的过程记录并同步引用。
仅当报告包含独立数值证据、未解决事项或必要历史背景时继续保留；写明日期和源码基准。
本机日志、大型输入/输出和环境不纳入文档提交。
