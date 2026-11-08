# RMCDHF 单节点 MPI：全程 D 状态诊断与优化候选

日期：2026-11-09。源码核查基线：`b8d6faf9c7097a50eb61096fa5934c27dfeec607`。
范围：`rmcdhf_mpi`、`rmcdhf_mem_mpi` 的单服务器大规模 CSF 计算。
状态：**初步源码分析完成；服务器运行根因未确定；采样设施、性能优化和 OpenMP 改造均未实施**。
尚未核实服务器可执行文件是否对应此源码基线。本文不构成实测性能或物理验收报告。

## 1. 现象与已知条件

以下来自用户反馈，尚无同步的服务器采样记录：

| 项目 | 已知情况 / 尚缺信息 |
| --- | --- |
| 程序 | 普通 MPI 和 memory MPI 均出现问题 |
| 规模 | CSF 超过 200 万时经常发生；总 CSF 数还是最大对称性 block 大小待核实 |
| 运行方式 | 通常单服务器；使用全部 48 个 CPU 核，实际 MPI rank 数、每 rank 线程数和绑定方式待记录 |
| 内存 | 512G 以上，用户确认容量足够；尚无作业 RSS、MemAvailable、swap 和作业限额的时间序列 |
| 存储 | 数据及 `MPI_TMP` 位于 NVMe；实际打开的文件和挂载信息待采样关联 |
| 表现 | 少数核心在运行，其他进程在 top 中呈 D；交替出现，可持续数分钟甚至数小时 |
| 轨道设置 | 内层轨道冻结，只优化新增外层轨道，例如 9s |
| 发生阶段 | **整个计算过程均会出现。** 曾在 `Lagrange multipliers:` 输出附近观察到，但用户没有持续监视，不能据此定位到 SETLAG |

最初把 Davidson 重复工作区和内存压力列为优先假设；结合后续硬件信息，
**不能继续把内存不足当作既定原因**。重复分配是源码事实，是否产生回收、换页或等待须由现场证据判断。
同样，NVMe 不等于已经排除所有内核 I/O 等待；设备类型本身也不能证明磁盘就是瓶颈。

Linux 的 D 表示不可中断等待，常见于 I/O 路径；仅看此状态不能区分文件读取、
回写、页等待或其他内核等待。普通 MPI 等待、负载不均、内存带宽饱和或线程超订阅，
本身都不足以解释持续的大量 D。应同时采集 D 与正在运行的任务，关联其内核等待位置和程序阶段。
依据：[proc_pid_stat](https://man7.org/linux/man-pages/man5/proc_pid_stat.5.html)、
[wchan](https://man7.org/linux/man-pages/man5/proc_pid_wchan.5.html)、
[内核栈](https://man7.org/linux/man-pages/man5/proc_pid_stack.5.html)。

## 2. 源码确认的工作量与并行边界

### 2.1 诊断必须覆盖完整计算链

[SCFmpi](../../src/appl/rmcdhf_mpi/scfmpi.f90) 在 EOL 路径先执行 MATRIX/NEWCO，
每轮再执行 SETLAG、按顺序更新未冻结轨道、输出轨道，随后重新构造矩阵和更新 CI。
[IMPROVmpi](../../src/appl/rmcdhf_mpi/improvmpi.f90) 包含系数/势构造、系数收集及径向求解。
因此应覆盖加载 MCP、矩阵构造、初猜/对角化、SETLAG、轨道更新和输出，
不能把最后一行 stdout 当作所有 rank 的当前执行位置，也不能仅围绕一次观察加计时。

### 2.2 SETLAG 的冻结侧势构造是独立优化候选

[普通 SETLAGmpi](../../src/appl/rmcdhf_mpi/setlagmpi.f90) 与
[memory SETLAGmpi](../../src/appl/rmcdhf_mem_mpi/setlagmpi.f90) 都先为轨道对 J、M
执行 `SETCOF → YPOT → XPOT → DACON`，保存两侧势，再根据 `LFIX` 选择公式。

- 当 J 冻结、M 可变时，最终分支只使用 M 侧势；J 侧势已经构造，但不参与该分支的最终表达式。
- 当 M 冻结、J 可变时，对应只使用 J 侧势。
- 两侧均可变时需要两侧势。普通版与 memory 版此分支的公式存在差异，必须各自保持，不能借优化统一算法。
- 对连续的 `1s–9s、2s–9s、3s–9s…`，现有 `MLAST` 已缓存相邻重复的 9s。
  这里的候选浪费是每次更换的冻结 J 侧构造，不能描述成每行都重新构造 9s。

可评估按分支延迟构造所需侧，进一步评估单次 SETLAG 内按轨道缓存。
实施前须审计全局工作数组、后续调用和调试输出的副作用；轨道或 CI 系数变化后缓存必须失效。
冻结轨道的波函数仍参与积分及正交约束，不能删掉其物理贡献。

每对乘子先完成一个标量 `MPI_Allreduce`，再由 rank 0 输出数值。
因此该行只证明 rank 0 已返回这次归约，不等于所有进程已经离开调用，也不定位下一次等待。
**上述冗余有源码依据，但尚未测量耗时，更未证明是全程 D 状态的原因。**

### 2.3 普通版文件扫描与 memory 版内存扫描

[普通 SETCOF](../../src/appl/rmcdhf_mpi/setcof.f90) 会 rewind 并扫描 MCP 记录，
按轨道筛选贡献、读取 IROW/系数并调用 DSUBRS 收缩。
[memory SETCOF](../../src/appl/rmcdhf_mem_mpi/setcof.f90) 仍有重复标签扫描、IROW 分配、
系数复制和收缩；[read1_mem](../../src/appl/rmcdhf_mem_mpi/read1_mem.f90) 与
[read3_mem](../../src/appl/rmcdhf_mem_mpi/read3_mem.f90) 从内存数组取数，名称中的 read 不代表这里执行磁盘读。

memory 版消除了这部分运行期文件读取，不等于消除了全程序 I/O、内存复制、分配或页等待。
两种程序都出现 D，不能直接归因于普通版 SETCOF 的 MCP 读取。
后续可比较两版在相同科学输入下各阶段的扫描量、复制量与耗时；数值结果仍以各自基线验收。

### 2.4 Davidson 与重复工作

| 源码事实 | 含义及待验证边界 |
| --- | --- |
| [MANEIGmpi](../../src/appl/rmcdhf_mpi/maneigmpi.f90) 每 rank 分配 `LWORK = 2*NCF*LIM + 2*LIM*LIM + 11*LIM + NVECT`，`LIM=min(NCF,2*NVECT+80)` | 工作区不随 rank 数自动均分；NCF 是当前 block 的规模，NVECT 是该 block 所选根的最高编号。不能把总 CSF 数和根数量直接代入 |
| [MATRIXmpi](../../src/appl/rmcdhf_mpi/matrixmpi.f90) 的 EMT/IROW 按本 rank 的 NELMNT 分配 | 不能说每个 rank 都保存完整稀疏 Hamiltonian；向量/工作区复制和矩阵分片是两回事 |
| [SPICMVmpi](../../src/lib/mpi90/spicmvmpi.f90) 按列轮转分工，但每次清零完整结果并归约 `N*M` 个 double | N=200 万、M=1 时应用缓冲区为 16 MB；这不是实际互连传输总量。列数相近也不保证非零元/计算量相近 |
| [gdsummpi](../../src/lib/mpi90/mpiu.f90) 分配临时 y、初始化、Allreduce 后复制回 x | 可测量临时缓冲与拷贝开销；原位归约等改写需独立回归，不能预先承诺收益 |
| [INIESTmpi](../../src/lib/mpi90/iniestmpi.f90) 合并初猜矩阵后各 rank 分别调用 DSPEVX；调用方上限为 4000 | 有重复求解候选，可评估单 rank 求解后广播，保留状态/错误传播；只解释相关初猜阶段，不能外推全程 |

量级示例：若**单个 block** 恰为 200 万 CSF、NVECT=1、double 为 8 字节，
上述 WORK 约为每 rank 2.44 GiB；假设启动 48 ranks，仅此项约 117.3 GiB。
这是分配公式估算，未包含 MCP、矩阵及其他数组，也不是实测驻留内存或内存不足结论。

## 3. MPI 优化与 OpenMP 的判断

当前优先顺序是：全程取证、消除已确认且通过回归的冗余、按阶段测量扩展性，再决定是否引入线程。
单节点 MPI 可以使用共享内存通信，但应用自己的数组不会因此自动变为共享。
Open MPI 的具体 transport 依版本与构建配置而异，不能直接套用某版参数；参见
[共享内存通信](https://docs.open-mpi.org/en/v5.0.x/tuning-apps/networking/shared-memory.html)。

先记录实际 BLAS/LAPACK 和线程设置，核实是否存在 48 ranks 叠加数学库线程。
线程总数本身不能证明这些线程同时活跃。若有多个 NUMA 域，再评估绑定和内存放置，
不要仅根据 48 核推定插槽拓扑。参见
[OpenBLAS 运行变量](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/) 与
[Open MPI affinity](https://docs.open-mpi.org/en/v5.0.x/tuning-apps/affinity.html)。

直接把全程序改成 OpenMP 需要重做状态所有权和并发写入，不能只加编译开关：
[TA/TB 工作区](../../src/lib/libmod/tatb_C.f90) 等模块状态共享，
[DMERGE](../../src/lib/lib9290/dmerge.f90) 会累加 `DC(IDY(I))`，直接并行外层列循环存在写冲突。
轨道更新顺序、归约顺序和 MPI 调用线程约束也须保留或明确重新验收。
未来可选择局部线程内核、减少 rank 并共享 Davidson 工作区等方案；
只有实测表明相应阶段占主导、设计解决线程安全后才开展。
设置 `OMP_NUM_THREADS` 不会自动把现有求解器变为混合并行版本，也不能保证消除 D。

## 4. 待实施的诊断与对照

### D0：无人值守的全程外部采样

目的：捕获同一时刻的程序阶段、运行 rank 与 D 状态任务，回答“在等什么、由谁触发、出现在哪些阶段”。
由于问题间歇出现且用户不持续监视，应先提供自动记录设施，而非依赖人工看到停顿后执行命令。
**本次只制定需求，尚无新增采样脚本。**

采样器应具备以下边界：

1. 按本次作业的启动器、进程树或 cgroup 限定对象，记录 PID/TID 与 rank 映射及进程起始时间，
   避免误采同用户的其他任务或 PID 重用；MPI 映射无法取得时显式标记未知。
2. 周期记录状态、CPU、驻留内存/缺页、I/O 计数、线程和绑定；同时记录设备延迟/队列、
   MemAvailable、swap 活动及 CPU/memory/io PSI。1–5 秒可作初始可配置周期，实际开销须验证。
3. 检测到 D 后对对应 `/proc/PID/task/TID/wchan` 和 `stack` 作限频采样，
   同时保留正在运行的 rank 样本。内核配置或权限不足时记录缺失，不自动提权。
4. 记录墙钟和单调时间，关联应用日志；仅凭外部日志不足以确定阶段时保持“未知”。
   采集实际输入/临时文件的路径和挂载关联用于现场诊断；原始机器路径和运行输出不提交仓库。
5. 日志有容量上限，任务结束后清理采样子进程；不修改计算输入，不主动清缓存或改变内核/MPI 参数。
   先做启用/关闭采样的对照，核查采样开销及是否扰动现象。

实现可使用 `/proc`、`ps -L`、pidstat、vmstat 和 iostat；缺少工具须降级并记录。
[pidstat](https://man7.org/linux/man-pages/man1/pidstat.1.html) 和
[iostat](https://man7.org/linux/man-pages/man1/iostat.1.html) 的时间序列应与等待栈一起解读，
不能只用 NVMe 的 `%util` 判断是否饱和。
[PSI](https://docs.kernel.org/accounting/psi.html) 用于关联资源等待；主机容量不能替代
[cgroup 作业限额](https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html) 的检查。

### D1：可关闭的每 rank 阶段计时

在 D0 不能充分定位时补充加载、MATRIX、INIEST/DVDSON、SETLAG、轨道更新及输出的阶段标记。
热点内部再按需拆出 SETCOF、YPOT/XPOT/DACON、局部矩阵乘和 collective。
记录调用次数、累计/最大耗时及实际工作量，在自然阶段边界汇总 rank 分布。
collective 内的时间包含等待慢 rank 的可能性，不能全记成通信传输成本。
不为测量额外插入 barrier，不在系数内循环逐条写日志；采用低频缓冲记录并验收诊断开销。

### D2：固定输入的单节点扩展性对照

先选能复现问题的输入与对应较小回归夹具，分别建立普通版、memory 版基线。
记录源提交/二进制 hash、GNU Fortran 15 具体版本、MPI/BLAS、编译配置、硬件拓扑、
每 block CSF 数、目标根、冻结集合、轨道顺序、迭代设置及科学开关。
在核实实际线程控制后，可比较 12/24/48 ranks 或先做 24/48 ranks；该矩阵尚未执行。

**每个 rank 配置须重新生成或使用经验证的工具重新分片 MCP，并核对内容。**
[SETMCP](../../src/appl/rmcdhf_mpi/setmcp.f90) 检查文件头的 MYID/NPROCS，
不能直接复用 48-rank 文件测试 24 ranks。
固定输入与初始轨道，记录初次读取和缓存已热的运行条件，分别报告，避免靠清系统缓存制造对照。

检查总墙钟、分阶段耗时、rank 工作量/等待分布、D 持续时间与等待栈、设备 I/O、内存和 PSI。
待检验的假设是：较少 ranks 或针对性的优化可能降低资源争用/冗余，改善总耗时。
结果可能不支持该假设；不能预设“满核最慢”或承诺加速倍数。
跨 rank 的归约舍入变化按事先登记的数值门槛验收，不能看到差异后放宽阈值。

### D3：按证据选择优化，分别回归

| 现场证据 / 独立源码候选 | 后续行动 |
| --- | --- |
| D 等待栈与设备读取、回写或页等待同期出现 | 定位具体文件/内核路径及诱因，比较普通/memory 阶段行为；不先重写并行模型 |
| 工作量分布明显偏斜，慢 rank 与 collective 等待相关 | 按非零元/贡献数分析分工；先区分计算偏斜与 I/O/调度造成的慢 rank |
| 活跃数学库线程超订阅、绑定或 NUMA 放置问题得到证实 | 固定线程与绑定做单变量对照；比较科学输出和墙钟 |
| SETLAG 冻结侧冗余 | 完成副作用审计后按需构造势；此优化独立于 D 根因，分别报告收益 |
| Davidson 复制/归约或 INIEST 重复求解占主要耗时 | 评估缓冲复用、归约改写或初猜单 rank 求解；涉及共享库时扩展回归范围 |
| 经上述优化仍由可并行内核或工作区复制主导 | 再设计局部 OpenMP / MPI+OpenMP，明确线程私有状态、冲突写入及 MPI 线程级别 |

SETLAG 至少覆盖冻结左侧、冻结右侧、两侧可变和重复轨道对，保留两版原有不同公式；
覆盖适用的 EOL/EAL 路径、调试行为、1/2/4 ranks 及完整一轮 SCF。
保持运算顺序的改写按同 rank、同工具链基线要求精确一致；改变归约顺序的方案预先登记容差。
检查能量、目标态/本征残差、乘子、径向输出、正交性和停止状态，遵循
[F2018 数值验收协议](f2018_migration_gfortran15_2026-11-08.md)。
普通版和 memory 版分别与自身基线比较，不能以两者数值强行相等作为验收。

## 5. 实施状态与后续交付

| 项目 | 2026-11-09 状态 |
| --- | --- |
| 用户运行条件和“全程发生”的修正 | 已记录，未以现场数据独立核实 |
| 源码调用链及冗余候选核查 | 已完成初步分析，未完成优化前副作用审计 |
| D0 自动采样、D1 阶段计时 | 待实施 |
| D2 服务器对照与 D 根因定位 | 未执行，未得到根因结论 |
| D3 性能修改及 OpenMP 方案 | 均未实施，未获得加速或物理验证结果 |
| 本次文档检查 | 核查源码引用与相对链接；不运行构建或科学计算 |

下一项交付应为 D0 采样设施及使用说明，然后用一轮可复现的服务器记录决定是否进入 D1/D2。
运行记录保留在本地或外部数据目录；完成后将源版本、输入 hash、采样覆盖、缺失项、
实际结果及结论边界写入 `docs/reviews/`，并更新本表。
此前的[严格构建与 NP 验证](../reviews/f2018_strict_interfaces_npol_2026-11-09.md)
仍按原范围有效，不能作为本问题已解决或已完成大规模性能验收的证据。
