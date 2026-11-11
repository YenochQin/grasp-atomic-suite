# 普通 rangular 的 OpenMP CSF 对计算

日期：2026-11-10。基于 `e568a45` 的工作树增量。
状态：首版 `rangular_omp` 已实施；串行 `rangular` 保留。Check/Release 小型回归和
23,459-CSF 单机对照通过，未完成全部态空间、所有平台或完整 SCF/RCI 物理验收。
当前生产 flags 已精简，CSFG 已隔离；现行配置以[构建策略](build_profiles.md)为准。

基于 `c5f36d2` 的 [2026-11-11 本机性能对照](../reviews/rangular_mpi_openmp_benchmark_2026-11-11.md)
补充 14,551-CSF 样本的三轮 MPI/OpenMP 及文本/描述符测量：33 次运行全部通过输出验证，
8 并发整程中位数为 MPI 4.06 秒、OpenMP 10.56 秒，适用范围和转换成本见该报告。

## 1. 入口与范围

`rangular_omp` 并行生成普通 Coulomb 角系数，支持原文本输入、描述符输入、完整相互作用及一阶限制。
输出仍为普通串行 `mcp.30…`，文件头 `MYID=0,NPROCS=1`。
线程数不改变 CSF 次序、稀疏位置、贡献次序、截断值、F0 过滤或积分标签约定。

排序、稀疏结构整理及最终文件写入仍串行。没有替换 MPI 目标，也没有改造 RMCDHF/RCI/CSFG。
多 rank 求解器所需的 MCP 分片仍须单独生成；不能把本目标的文件改名后直接交给多 rank RMCDHF。

构建开关 `GRASP_RANGULAR_OPENMP` 默认跟随初次配置时的 `GRASP_JJ2LSJ_OPENMP`，也可独立指定。
`configure.sh --no-openmp` 同时关闭两项；启用 angular OpenMP 必须找到 Fortran OpenMP，不能静默退回串行。
关闭该开关时不构建 `rangular_omp`，普通 `rangular` 不依赖 OpenMP runtime。

## 2. 构建与运行

从 suite 根目录执行；构建目录应尚不存在：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC $FC

./configure.sh --serial --no-mpi --no-csfg \
    --profile f2018-release --build-dir build-rangular-omp \
    -DGRASP_JJ2LSJ_OPENMP=OFF -DGRASP_RANGULAR_OPENMP=ON

cmake --build build-rangular-omp \
    --target rangular rangular_omp csf_descriptor_verify -j4
```

在独立计算目录运行，仍按原流程输入完整/一阶相互作用选项：

```fish
set -gx OPENBLAS_NUM_THREADS 1
/path/to/suite/build-rangular-omp/bin/rangular_omp \
    --threads 8 --csf-descriptors /path/to/csfs.bin \
    --timings timings.csv
```

不指定 `--csf-descriptors` 时默认读当前目录 `rcsf.inp`，也可使用 `--csf-text FILE`。
原有路径校验与新建计时文件规则见[输入说明](rangular_descriptor_input.md)。

| 选项 | 默认与含义 |
| --- | --- |
| `--threads N` | 默认使用 OpenMP 环境设置的最大线程数并受 thread limit 限制；显式值必须为正且不超过 limit。 |
| `--pair-chunk N` | 每个任务最多 256 个连续 CSF 对；行末会缩小任务，使可用任务分散给线程。 |
| `--pair-buffer-mib N` | 所有待提交二电子结果的总缓冲预算，默认 64 MiB；0 强制走独立临时流，主要用于回归。 |

缓冲预算不包括共享 CSF 数组、每线程当前对工作区、任务元数据、原 block 系数累积和排序数组。
它不是进程 RSS 上限。普通 `rangular` 拒绝大于 1 的 `--threads` 或并行缓冲定制，避免误以为已启用并行。

用[版本对照运行器](rangular_version_comparison.md)时，显式指定候选线程数：

```fish
../graspkit-tools/.venv/bin/python scripts/compare_rangular_versions.py \
    --reference /path/to/original/rangular \
    --candidate /path/to/suite/build-rangular-omp/bin/rangular_omp \
    --candidate-threads 8 \
    --csf /path/to/original.c --binary /path/to/csfs.bin \
    --verifier /path/to/suite/build-rangular-omp/bin/csf_descriptor_verify \
    --output /path/to/new-results --repeat 1
```

运行器将原版及数学库固定为单线程；候选的 `OMP_NUM_THREADS` 和 `--threads` 使用指定值，
记录在 summary/execution JSON 中。不指定 `--candidate-threads` 时，运行器仍按单线程基线执行。

## 3. 状态与构建边界

并发闭包从 `ONESCALAR`、`RKCO_GG` 和实际 `CORD` 回调展开，共 97 个过程，
清单为 [concurrent_sources.cmake](../../src/appl/rangular90/openmp/concurrent_sources.cmake)。
[机械审计](../../test/rangular/audit_thread_closure.py)检查调用覆盖、模块依赖和实际编译命令；它不是形式化线程安全证明。

| 数据 | 所有权 |
| --- | --- |
| 当前 block 的 IQA/JQSA/JCUPA、轨道、控制参数及初始化后的只读表 | 共享只读；加载和释放均在并行区域之外。 |
| `m_C`、`dumx_C`、`trk_C`、`buffer_C` | OpenMP 定向对象中的 threadprivate 数据，每个 worker 独立使用。 |
| `SIXJ` 的 `RACA` | 现有 threadprivate lazy cache，在并发对象中启用。 |
| `memory_man` 的局部临时变量及并发闭包局部数组 | 并发对象用自动存储；不复用按 `-fno-automatic` 编译的分配器/计算对象。 |
| 结果块 | 每个任务独占写，barrier 后由初始线程读取；只复制数据，不跨线程传递 TLS 指针。 |
| `sacoef_C`、全局 npos、临时结构文件和最终 MCP | 初始线程独占提交，沿用原串行整理逻辑。 |

四个共享状态模块仅增加 OpenMP sentinels。普通库将它们视为注释；新目标从原源文件定向重编译
四个状态模块、内存管理模块和并发调用闭包。直接对象优先于旧静态 archive，
同一可执行程序只链接一份模块定义，没有复制整套 libmod/lib9290/librang90。
Linux ELF 检查确认 `m_C/dumx_C/trk_C/buffer_C` 与 RACA 在新目标中为 TLS。

并发对象与调用端的 `.mod` 目录独立，threadprivate 接口优先于原目录。
整个并发闭包使用 `-fopenmp`，由 GNU OpenMP 隐含自动/可重入局部存储。
全局隐式 SAVE 配置已删除，不再添加 `-frecursive` 或 `-fautomatic` 兜底。
RACA 的持久缓存语义用源码 SAVE 明确表达，在 OpenMP 对象中保持线程私有。

每个 block 建立固定大小的线程组，并关闭动态改变线程数。工作区在 block 开始时分配，结束时释放；
各批次使用相同线程数和 affinity policy，无嵌套 parallel、可迁移 task、teams 或 resource pause。
这些约束满足 [OpenMP threadprivate 持久性条件](https://www.openmp.org/spec-html/5.2/openmpse24.html)。
该适配不意味着整个 rang90 已成为可嵌套调用的通用 context 库。

## 4. 任务、内存及输出

[并行控制器](../../src/appl/rangular90/openmp/rangular_parallel.f90)按当前 JA 列切分 JB 范围，
每批最多四倍线程数个任务，使用 `schedule(dynamic,1)` 分配。选定一阶限制时，先处理参考区间，再处理对角项。
任务数和结果缓冲不随全部三角区域预分配，不建立 NCF² 大小的任务表。

每个任务在同一 worker 中按原顺序执行该对的 ONESCALAR、RKCO_GG，保存一电子结果和原始二电子贡献。
ONESCALAR 早退时只依赖其已定义的第一轨道标签；拒绝的对不读取未定义的第二标签或系数。
一个结果块满后，将其已有和后续贡献转写独立 `newunit` scratch stream；用 int64 字节位置读回，正常结束时删除。
只复制 LABEL 的五个有效 Coulomb 字段，未使用的第六槽不进入结果队列。

barrier 后，初始线程按原 JA/JB 次序取得结果，交回原 [MCP](../../src/appl/rangular90/mcp_gg.f90) 处理。
因此位置编号、重复贡献、浮点计算顺序及原排序输入顺序都保持，工作领取顺序不决定文件顺序。
没有跨线程浮点 reduction，也没有用一个大 critical 区锁住整个角系数计算。

新的缓冲/临时流错误在批次边界汇总并非零退出。深层旧库仍保留原有 fatal STOP 行为；
使用 fresh 输出目录及对照运行器的文件完整性检查判定成功，不只看旧例程的退出码。
显式角动量 debug 输出不在首版并行范围，检测到时拒绝运行并提示使用串行目标。
GNU runtime 的并发调用边界见[线程安全说明](https://gcc.gnu.org/onlinedocs/gcc-15.2.0/gfortran/Thread-safety-of-the-runtime-library.html)。

## 5. 验证与实测

[七项线程回归](../../test/rangular/test_openmp.py)覆盖：

- 1/2/4/8 线程、真实多 worker 工作量、多 block 工作区重建，完整和一阶限制。
- 文本/描述符输入、不同 chunk 大小、强制结果流回退。
- 四个活跃壳层、f 壳 seniority、较高 j 的表查询及原 jj2lsj 测试基组。
- CLI 正常/非法参数、候选线程数传递、数学库保持单线程。

每个组合比较完整 MCP 文件 SHA-256，要求逐字节一致。
按最新范围和 flags 收敛后，普通 Check/Audit/生产 Release 各 44 项通过；关闭 OpenMP 的 Check
42 项通过，普通 MPI 生产配置 78 项通过。实际编译清单中没有 CSFG 源文件。
下面的 45/43/98 数量是范围收敛之前的阶段记录，不作为当前 CSFG 验收结论。
独立 angular OpenMP 的严格 Check/Release 各 45 项 CTest 通过；禁用 OpenMP 的 Check 43 项通过。
MPI/CSFG/jj2lsj OpenMP 全开 Release 98 项通过；后续输入兼容修复又通过相关专项回归。
旧 audit 存储配置也曾通过专项测试；按后续范围收敛决定，该兼容路径现已移除，不作为当前支持目标。

audit 复核曾捕获 `--help` 字符串从函数返回后长度变零，以及动态 basename 赋值失败。
现采用显式字符串输出参数；后续直接移除了隐式静态局部存储及其抵消选项。
临时诊断已移除，保留正向 help 与多 block 输入回归。

大样本目的：在固定输入和 O2 数学设置下，验证线程数变化保持全部 MCP 内容，同时测量单机扩展性。
输入为用户提供的 23,459 CSF、5 block（最大 6,276）、41 轨道；候选对 63,065,280，
每次输出 38,461,197 项角系数。使用同一二进制 CSF、完整相互作用、数学库单线程，当前环境可见 8 个 CPU。
预期是字节一致，不预设加速阈值。每种设置仅运行一次，文件缓存未清空。

| 执行方式 | 整程墙钟 | 相对串行 | 相对 OpenMP 1 线程 |
| --- | ---: | ---: | ---: |
| 原串行控制路径 | 67.64 s | 1.00× | — |
| OpenMP 1 线程 | 64.51 s | 1.05× | 1.00× |
| OpenMP 2 线程 | 41.78 s | 1.62× | 1.54× |
| OpenMP 4 线程 | 30.65 s | 2.21× | 2.10× |
| OpenMP 8 线程 | 24.50 s | 2.76× | 2.63× |

五次运行 MCP 均逐字节一致。8 线程的成对生成与缓冲约 14.72 s，其中并行批次约 11.38 s；
串行排序/输出约 8.72 s，结构整理约 1.02 s，限制了整程扩展。
1 线程与原串行还有调用组织/编译差异，因此另外给出相对 OpenMP 1 线程的扩展性。
本表为精简 flags 前的首轮实测；精简后的三轮 8 线程对照见[构建策略](build_profiles.md)。
本机结果不能直接外推其他 CPU、NUMA、线程数或更大基组。

可追溯记录位于忽略的 `build-rangular-omp-release/audit/scaling-23459-v1/`，包含 summary、
输入/可执行文件 SHA-256、日志、逐阶段 CSV、全部 MCP 和实际测量的可执行文件副本。
测量后继续修复输入生命周期、简化启动诊断并整理 flags；角系数内核、任务调度和输出算法保持。
本轮没有进行完整 SCF/RCI、活性空间收敛或物理精度验证。

## 6. 剩余边界与警告审计

新控制器没有编译器诊断；定向编译仍会再次报告旧源中的诊断。
Check 闭包含 58 项原有 real→integer 转换诊断，主要涉及半整数/相位编码；表达式保持。
Release 另外有 RMEAJJ 的小-j term 编码前置条件、NINE0 的零参数前置条件及 SETSDA 非空结构前置条件诊断。
SETSDA 在此路径始终保留所有对角项；小-j term 编码来自原表。NINE0 的旧零参数判定及极端整数范围仍属于遗留数值审计范围。
没有把新增目标的重复诊断计成新修复，也没有宣称遗留正确性告警全部清零。

警告预算更新记录在 [baseline JSON](../reference/f2018/gfortran15-warning-baseline.json) 的
`rangular_openmp_budget_update`。独立 angular 开关用 `-a1/-a0` 区分非默认组合；
旧四轴含 CSFG 的矩阵停止作为本轮工作范围。2026-11-11 基于 `3f902c2` 补齐独立开关后，
当前普通矩阵包含三种 profile，各 12 个有效组合（共 36 个）；serial 关闭时不重复枚举无目标的 angular 开关。
生产 profile 不进行警告预算认证；24 条预算覆盖全部普通 Check/Audit 组合，
包括只开启 jj2lsj OpenMP、只开启 rangular OpenMP 的 MPI/非 MPI 配置。
回归测试同时检查独立开关枚举和已提交预算覆盖，避免可构建的配置在验收时找不到基线。
CSFG 与旧 Release 诊断预算保留在 Git 历史中，**没有据此宣称完整 36 组合矩阵验收**。

当前范围不含跨 block 并行、并行排序、MPI 输出重分片或全库显式 context 重构。
任务缓冲预算不替代旧 MCP 整数容量和 block 累积数组限制；高并发、内存耗尽及所有罕见态还需各自验收。

## 7. 独立开关修复的本机 CSF 对照（2026-11-11）

基于 `3f902c2` 的验收矩阵/预算增量，在 macOS arm64、GNU 15.3.0 上使用用户提供的
14,551-CSF 文本样本，四个 block 分别为 1,212、3,298、4,765、5,276 个 CSF。
目的为验证独立开关构建的程序在不同线程数下保持完整 MCP 输出；预期为逐字节一致。

使用 Check 配置、完整相互作用、数学库单线程。串行控制程序来自 jj2lsj OpenMP 开启、
rangular OpenMP 关闭的配置；候选来自相反的开关组合。两个配置均启用普通 MPI 构建，
此次运行的是普通串行/线程程序。输入和程序均冻结为独立副本，重新生成串行基准，
按 1/2/4/8 线程运行候选，比较零容差数值结果及完整文件 SHA-256。

四种线程设置均正常结束，10 个 MCP 文件全部与串行基准逐字节一致，符合预期。
证据保存在忽略目录 `build-review-openmp-independent/zn-csf-frozen/`：
`summary.json` 记录输入/程序身份、线程设置和每次运行，`comparison-threads-*.json`
记录全部 MCP 文件哈希，各运行目录保留输入、输出、日志及阶段计时。
运行期间有其他编译/测试活动，因此耗时不作为扩展性结论。
本次未验证完整 SCF/RCI 或物理收敛；完整 MPI 验收仍需在 Linux 补跑
[构建策略](build_profiles.md)所列的四项平台专属测试。
