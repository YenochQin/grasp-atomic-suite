# rangular：描述符输入与单节点 OpenMP 的可行性及实施边界

日期：2026-11-10。源码基准：suite `e229657d8a32b851da9c130928545bcb3ac209f1`；
rCSFs `6e4fcff1439fb1b56d4891908c5bd58722ddf64a`。

现行状态：**描述符输入适配和计时已实施；普通 rangular 的首版 OpenMP 现已实施并完成所述回归**，
见 [OpenMP 现行说明](../implemented/rangular_openmp.md)。采用审计后的 threadprivate 过渡适配；全库显式 context 尚未完成。
实现、命令、已通过的软件回归及支持范围见[描述符输入说明](../implemented/rangular_descriptor_input.md)。
下文保留最初的设计分析及后续方案，不能将其每项候选均视为已实施；最初研究仅新增本文和索引，当时已有的 MCP 工具记录不算作该次运行结果。

适用范围：普通 `rangular` / `rangular_mpi`，不包含 CSFG 专用表示。
本文承接 [rangular 数据流](../implemented/rangular_flow_and_state.md)、
[F2018 迁移计划](f2018_migration_gfortran15_2026-11-08.md)和
[rCSFs 接入讨论](../../../rCSFs/docs/grasp_descriptor_integration_discussion_20261110.md)。
最后一份讨论中希望避免 rCSFs 编译依赖的取向继续有效。

## 1. 判断与建议边界

两个方向在技术上都可行，但应独立实施、分别验收：

1. **描述符输入**：physical-seniority Int8 描述符与 GRASP 的 CSF 数据已经接近；增加确定性的输入适配层，
   可绕过“恢复文本、再解析文本”。主要收益候选是数据准备与加载，不能由此推断角系数生成也会同样加速。
2. **OpenMP**：从串行程序的 CSF 对计算出发适合单节点目标；工作重点是完整调用链的可变状态隔离及输出组织，
   不是把 MPI 调用逐个改成 OpenMP 指令。收益、线程数和等待原因必须实测。
3. 两者可以共用一个明确的边界：**共享只读的当前 block CSF 基组 → 独立角动量工作区 → 有序输出**。
   描述符能在原串行内核上验证，OpenMP 也能先继续使用文本输入；两项互不构成开始工作的前提。

现有 F2018 计划将新增计算线程排除在其原范围之外。这里作为新增性能工作线记录，
不把线程改造计入“语言标准迁移已完成”，也不必等待全套程序迁移结束才开始这个局部试点。
进入并发前，rangular 实际调用闭包中的接口、状态寿命和未定义行为必须单独完成检查。

## 2. 输入接口：需要适配什么

### 2.1 数据布局与语义

rCSFs 的权威契约是 [descriptor_schema.rs](../../../rCSFs/src/descriptor_schema.rs)，
编码/恢复实现在 [descriptor_v2.rs](../../../rCSFs/src/descriptor_v2.rs)，
seniority 表在 [seniority.rs](../../../rCSFs/src/seniority.rs)。
suite 的真实落点是 [LODCSH2](../../src/lib/lib9290/lodcsh2.f90)，
并由 [IQ](../../src/lib/lib9290/iq.f90)、[JQS](../../src/lib/lib9290/jqs.f90)、
[JCUP](../../src/lib/lib9290/jcup.f90)、[ITJPO](../../src/lib/lib9290/itjpo.f90)、
[ISPAR](../../src/lib/lib9290/ispar.f90) 消费。

| 信息 | 描述符侧 | suite 所需适配 |
| --- | --- | --- |
| 占据数 | `sub_n` | 对齐 header 的轨道次序，填 `IQA(i,icsf)`；补齐固定 core。 |
| seniority | `sub_v` 为物理 `v` | 填 `JQSA(i,1,icsf)`，无需恢复旧打印标签再查表。 |
| 额外态标签 | 当前四通道没有 `w` | `JQSA(i,2,icsf)` 是 `w`。第一阶段仅支持确认无需额外标签的态；不猜测简并标签。 |
| 子壳层角动量 | `sub_2j`；`-1` 表示未打印 | 补全物理值后存 `2J+1` 到 `JQSA(i,3,icsf)`。空/满壳层为 `1`；不能把 `-1+1=0` 当答案。 |
| 中间耦合 | `sub_2k`，与描述符子壳层槽位关联 | 补全省略耦合，按开放壳层链压缩到 `JCUPA(1:nopen-1,icsf)`，值为 `2K+1`；不是原数组整列加一。 |
| 总 J 与宇称 | 尾部 `total_two_j`、`parity` | 在 `JCUPA(NNNW,icsf)` 存 `(2J+1)*parity`；这是特殊保留槽。 |
| 基组身份 | header、轨道表、行序、对称性信息 | 保留 core/peel、`NP/NAK/NKL/NKJ/NH`、电子数、block 边界、block 内 CSF 顺序及原行映射。 |

上述赋值可定位到 `LODCSH2` 的 186–190、281–317、340–347 行。
`J=0` 的部分占据壳层仍是开放壳层，不能因角动量为零就从开放壳层计数中删除。
从全轨道描述符到开放壳层链，需要区别空壳、满壳、部分占据且 J=0，以及确实省略的中间耦合。

当前 [PACK](../../src/lib/lib9290/pack.f90) 实现直接赋值 byte 数组；旧注释不能作为仍需 bit packing 的依据。
[orb_C](../../src/lib/libmod/orb_C.f90) 和 [stat_C](../../src/lib/libmod/stat_C.f90) 中这三组运行期数组本来就是 byte，
因此不能把 Int8 输入描述成“将现有 GRASP CSF 内存从 Int32 降低四倍”。
计数、偏移、尺寸乘法和稀疏位置仍需宽整数及边界检查；描述符 `2J<=126` 只保护量子数表示。

### 2.2 支持范围不能由“v 已经一致”推出

[terms_C](../../src/lib/libmod/terms_C.f90) 的表项是 `(v,w,2J+1)`。
例如 `j=9/2,q=4` 行存在相同 `(v,J)`、不同 `w` 的态，而当前描述符没有独立 `w` 通道。
普通 `LODCSH2` 当前还将第二槽恒写零，显式 seniority 文本分支主要针对 `j=7/2,q=4`。
这些事实意味着：**描述符能编码、文本能恢复、此普通 reader 能正确表达全部态，是三个不同命题**。

第一阶段应以实际输入和对照夹具定义支持集合。若态不能由现有字段唯一指定，或与普通 reader 的选择不一致，
应给出明确的 unsupported 错误；是否增加 `w`/state identifier 属于单独的格式与角动量内核能力变更。
不能用固定 `w=0` 声称所有 GRASP 态均已支持，也不能把旧 reader 的限制当作物理正确性 oracle。

### 2.3 推荐的轻量输入路径

```mermaid
flowchart LR
    P[已有描述符 Parquet 与 header] --> S[suite 预处理脚本：分批解码与校验]
    S --> B[带版本的二进制交换文件]
    B --> F[Fortran：语义补全与数组装配]
    T[原文本 rcsf.inp] --> L[原 CSLH / LODCSH2]
    F --> C[同一 CSF block 数据边界]
    L --> C
    C --> K[原串行角系数内核 / 后续 OpenMP 内核]
```

预处理脚本放在 suite，复用现有 Python 环境中的可用 Parquet 库，实施时再核对依赖。
它不导入 rCSFs、不链接 Rust 库；保留原始 `v/2J/2K` 及缺省标记，
由 suite 的 Fortran 输入层统一完成 GRASP 语义适配，避免在 Python 和 Fortran 中各写一套物理态补全。
描述符格式依然由 rCSFs 定义，通过固定兼容性夹具检测格式漂移，不在 suite 维护另一套生成算法。

Fortran 用 stream unformatted 输入读取约定字节。交换文件至少定义：magic、版本、端序、
字段宽度/符号、数组顺序、轨道及 core 表、每个 block 的长度/偏移/对称性、行序与来源标识。
计数和文件偏移采用明确的宽度；进入旧 default-integer 内核时检查范围并拒绝溢出。
需要检查 metadata 的 `seniority_encoding=physical`、`descriptor_dtype=int8`、`max_two_j=126`，
以及格式版本、header 关联、缺列/空值/非法宇称/截断和尺寸乘法溢出。

`access='stream'` 只消除顺序非格式化文件的记录标记，不自动定义跨平台字节序。
首版应明确支持的端序及转换/拒绝策略，而不依赖 GNU 扩展选项隐式处理。
大数据按批转换、逐 block 装配；不要额外保留整份 Parquet、全部解码数组和运行数组三套副本。
它是输入交换文件，不替代现有 MCP 格式。

接入必须同时提供 [CSLH](../../src/lib/lib9290/cslh.f90) 的 header/block 元数据入口与
`LODCSH2` 的逐 block 数据入口；仅替换后者，前者仍会扫描文本。
完整 header 中的 block lengths 必须保留，不能只按相同 J/宇称重新分组。
尤其一阶相互作用的 [GETINF](../../src/appl/rangular/getinf.f90) 使用前 `1…ICCUT` 个 CSF 定义参考区间，
所以输入适配不能擅自排序、去重或改变 MR 所在序号。

Parquet 有列块、row groups 和文件 metadata，不能按连续 Int8 文件直接读取；结构依据见
[Apache Parquet 官方格式](https://parquet.apache.org/docs/file-format/)。
以后如确有必要可增加原生导出或 C ABI，但不作为本阶段构建前提。

## 3. 当前串行/MPI 差异改变了性能比较方式

以下均来自当前实现，尚无新的计时结果：

| 事实 | 对设计/比较的含义 |
| --- | --- |
| 串行 [MCP](../../src/appl/rangular/mcp_gg.f90) 100–111 行分配 block 系数数组，168–211 行展示 T 的扩容及磁盘回退，349–367 行选择 SORTMEM/SORT | 串行已有内存优先路径；仍写稀疏结构 `tmp.30`，最终 MCP 也必须输出，不能称为完全无 I/O。 |
| MPI [MCPMPI](../../src/appl/rangular_mpi/mcpmpi_gg.f90) 149、197、206 行直接写 T、V 与结构临时记录，随后使用 SORT | MPI 与串行不仅分工不同，系数存储策略也不同。OpenMP 版若保留串行策略，收益须拆分为 I/O 与并行贡献。 |
| MPI 124–126 行按 `JA=1+rank,NCF,nprocs` 轮转列，列内遍历 `JB=1…JA` | 已不是连续等长列分段，不能仅凭三角区域认定必然严重失衡；选择定则、开放壳层结构和贡献数量仍会造成不同列代价差异。 |
| [GENMCPMPI](../../src/appl/rangular_mpi/genmcpmpi.f90) 148–152 行各 rank 分配三组 CSF 数组；[LODCSLMPI](../../src/lib/mpi90/lodcslmpi.f90) 由 rank 0 加载后广播 | OpenMP 可共享一份当前 block CSF 输入；这不会自动消除角系数输出、排序和线程局部工作区的内存。 |

单份三组输入数组按当前分配形状约为 `5*NNNW*NCF` 字节（byte 为 8 位时），不含其他工作区。
MPI 有 rank 级副本；OpenMP 的候选内存模型是：

```text
一份 CSF 基组 + 线程数 × 局部角动量工作区 + 有界待提交记录 + block 输出/排序工作区
```

不能给每个线程复制整个 `sacoef_C`，否则会抵消共享输入的意义。
现有 `NPOS/LLIST*` 和 MCP 文件的整数宽度还是规模限制；新输入文件支持 int64 计数不等于旧内核和 MCP 自动支持超大计数。

单节点 MPI 本身可使用共享内存通信；不是所有 MPI 消息都会走网络。
应用层数组复制和通信 transport 是不同问题，参见
[Open MPI 共享内存文档](https://docs.open-mpi.org/en/v5.0.x/tuning-apps/networking/shared-memory.html)。

## 4. OpenMP 的实施边界

### 4.1 先做可重入内核，再加入调度

推荐先把“一个 CSF 对生成有序贡献列表”抽成可独立调用的计算边界。
共享数据只包括当前 block 的基组、轨道信息、筛选设置和已初始化只读表；
每个执行上下文独立拥有随 CSF 对变化的工作区、缓存及系数缓冲。

已确认的隔离对象包括：

| 状态 | 源码证据及处理 |
| --- | --- |
| `m_C` | [SETQNA](../../src/lib/lib9290/setqna.f90) 写入 NQ/JJQ/JLIST/KLIST/JJC/NCORE/NPEEL；整个当前对工作区需独立。 |
| `dumx_C` | [ONESCALAR](../../src/lib/librang90/onescalar.f90) 保存/重排耦合链；不能只私有化循环变量。 |
| `trk_C` | [PERKO2](../../src/lib/librang90/perko2.f90)、ONESCALAR1/2、EL 分支使用 ID/IK/BD/BK 工作数组；也在并发闭包内。 |
| `buffer_C` | [SPEAK](../../src/lib/lib9290/speak.f90) 更新 NVCOEF、LABEL、COEFF 并可能调用 ALCBUF 扩容；指针和目标内存都需独立拥有。 |
| 6-j 缓存 | [SIXJ](../../src/lib/librang90/sixj.f90) 的 DATA 初始化 RACA 会按需写入；已有 threadprivate 指令仅在该文件以 OpenMP 编译时生效。 |
| 文件与日志 | SPEAK 等深层例程可能写单元 99，错误分支可能 STOP；应缓存到上下文或在首版明确采用串行调试路径，统一错误收集。 |

这不是完整线程安全清单。还需覆盖 `RKCO_GG → EL1…EL5` 的所有可达分支、
函数调用和回调，以及 DATA/SAVE/首次初始化标志；不能只从外层 USE/CALL 文本推出闭包安全。

首选显式 `angular_workspace` 参数及初始化/销毁过程，串行包装器保留旧接口，实现在原有共享库内。
可在小型过渡原型中使用完整审计后的 threadprivate，但其线程绑定、跨 parallel region 初始化、
指针分配与释放必须有规则，不能与可迁移任务随意混用。
`firstprivate(pointer)` 只复制关联，不能替代独立分配。
给 RKCO_GG 外面套大锁虽可避免同时进入，但会将主要计算再次串行化。

依据 [OpenMP 对被调用例程的数据共享规则](https://www.openmp.org/spec-html/5.2/openmpsu34.html)，
模块关联及 SAVE/DATA 状态不会因为外层使用 `default(none)` 就自动变为私有。

### 4.2 有界动态任务与确定性提交

初版保持逐 block 处理、串行排序；在 block 内按 JA 列或小段列领取工作。
每列内仍按原 `JB` 顺序计算，必要时将很重的列拆成固定 JB 片段。
候选调度从动态领取小任务开始，与轮转静态分工做实际对照；不预设某个 chunk 最优。

推荐数据流：

1. 调度器给每段定义确定的 `(JA,JB_begin,JB_end)` 及逻辑顺序，限制同时在途的记录总量。
2. 工作线程只计算本段，结果记录 `(JA,JB,通道,积分标签,贡献序号,数值)` 及结构存在标志。
   保留当前 CUTOFF、F0 过滤、轨道标签规范化和同一对内的贡献次序。
3. 按原 `(JA,JB)` 遍历顺序提交。使用每段非零对计数的前缀和或单一提交者分配 `npos`，
   对角元素始终保留，多个贡献共用同一稀疏位置；不能让线程直接竞争全局 `npos`。
4. 复用既有 block 累积及排序协议，再统一写 MCP。首版不同时并行化 SORTMEM/SETSDA。

可用同步批次先验证正确性，再视实测开销实现生产者/消费者流水线。
批次需要输出字节上限及背压，不能只限制任务个数后假定每列贡献数固定，
也不应为所有 `NCF*(NCF+1)/2` 对预分配记录。
一个很慢的早期任务会阻塞有序提交；这是需要测量的代价，细分任务可以缓解但有调度成本。
独立对计算不必引入跨线程浮点 reduction，因此有条件争取保持同一对的原运算顺序。

### 4.3 构建与下游文件接口

候选目标名可为 `rangular_omp`，以普通串行目标为基线；名称在实施前确定，不宣称已存在。
并发调用闭包必须统一采用可重入编译及正确 OpenMP 设置，不能只给 main 加 `-fopenmp` 后
继续调用按 `-fno-automatic` 编译的旧对象。
[GNU Fortran 15 文档](https://gcc.gnu.org/onlinedocs/gcc-15.2.0/gfortran/OpenMP.html)说明
`-fopenmp` 隐含 `-frecursive`，大局部数组还需注意线程栈；大的持有型工作区宜显式分配。

当前 [FortranProfiles](../../cmake/FortranProfiles.cmake) 的 legacy/audit 保留静态存储；
strict profiles 不保留，也仍不证明模块状态已可重入。
可参考 [jj2lsj CMake](../../src/appl/jj2lsj/CMakeLists.txt) 对并发闭包定向编译的组织方式，
但其闭包比 rangular 小，既有结果不能作为这里的验收证据。
共享库保持唯一源码所有权和全局状态定义，避免同一可执行程序出现两套模块实例；新增/关闭 OpenMP 都要构建回归。

**输出分片是一个独立接口决策。** 首版产生普通 `mcp.*` 最容易与串行基线比较，
但现有 [RMCDHF MPI SETMCP](../../src/appl/rmcdhf_mpi/setmcp.f90) 121–127 行要求文件头的 rank/进程数相符。
不能将 `mcp.*` 改名为 `mcp000.*` 就直接给多 rank RMCDHF 使用。
若后续仍运行 MPI 求解器，需要另行实现并验证“按目标 rank 数输出/重分片 MCP”，
重建轮转列归属、局部位置、结构与文件头；输出 rank 数与 OpenMP 线程数应独立。
若整条求解链也采用线程模式，则其改造是后续程序的单独工作，不由 rangular 试点自动完成。

## 5. CPU 等待需要单独取证

“单服务器 + CPU 空闲”支持尝试线程方案，但不是负载失衡或 MPI 通信过慢的证明。
已有 [RMCDHF 单节点诊断记录](rmcdhf_mpi_single_node_diagnostics_2026-11-09.md)
描述的是 RMCDHF 全程现象；不能直接把它定位到 rangular。
若观察到进程 D 状态，它指向不可中断等待，需查等待位置和 I/O/页活动，
而不是仅调度算法，见 [Linux proc_pid_stat](https://man7.org/linux/man-pages/man5/proc_pid_stat.5.html)。

对 rangular 建议记录：header/block 加载、ONESCALAR/RKCO_GG 成对计算、
扩容/临时写入、SETSDA、SORTMEM/SORT、最终写入的墙钟；另记录每 worker 任务数、
候选对/有效对数、T/V 贡献数、计算与提交等待时间、峰值内存、I/O 字节。
汇总在自然阶段边界完成，不在每个 CSF 对写日志或为计时增加 MPI barrier。
阶段总时间和嵌套子计时分开，不能重复相加。

输入优化的加速上限受加载占总时长的比例限制；如果读取仅占很小一部分，优化它对整程收益有限。
如果 I/O 或排序占主导，只并行角系数生成也会出现相应上限。
NUMA、绑定、数学库线程及物理核/SMT 情况由服务器实际配置决定，首轮避免并行层叠加造成超订阅。
不根据“内存容量大/NVMe”预先排除带宽或 I/O 等待，也不据源码分配量认定发生内存不足。

## 6. 建议的实施顺序和验收

| 阶段 | 交付 | 进入下一阶段的证据 |
| --- | --- | --- |
| A：基线与测量 | 固定输入/交互设置、工具链及二进制身份；串行和 MPI 分阶段记录 | 正确性基线可重复，计时开销已检查；服务器等待根因仍可保持未知。 |
| B：输入适配 | 交换格式、预处理、Fortran loader；先用原串行内核 | 同一输入的有效内部数组、header/block/顺序精确一致；非法/不支持态明确拒绝；生成 MCP 对照通过。 |
| C：状态重构 | 显式工作区及结果边界，仍单线程 | 与原程序一致；两个独立 context 交错运行、连续初始化/释放、失败恢复均通过。 |
| D：线程原型 | 有界动态调度及有序提交；先保留原输出/排序 | 同一工具链 1/2/4/8 线程、多次重复及不同调度下结构/系数一致；扩容、回退和错误分支覆盖。 |
| E：单节点评估 | 按实际核数与 NUMA 拓扑扩展线程/rank 对照 | 报告总墙钟、阶段占比、内存、I/O 和负载分布；依据结果再选默认线程数或优化排序/输出。 |

输入和线程特性保持正交开关，至少比较：文本串行、描述符串行、文本 OpenMP、描述符 OpenMP。
性能实验另记录一次性 Parquet 转换成本与缓存文件复用时的加载成本，缓存冷热条件分别报告。
对比 MPI 时补充各程序存储路径，避免把所有差异写成 MPI/OpenMP 的差异。

测试输入至少包含：空/满壳、单开放壳、部分占据 J=0、多个开放壳、显式零与缺省值、
`j=7/2,q=4` 不同 seniority、超支持态/歧义态、正负宇称、多 block、一阶限制及完整相互作用、
零贡献非对角对和仅结构对角元。文本 reader 的未使用 padding 可能未定义，
比较有效槽及保留槽，不把未初始化 padding 的随机字节作为 oracle。

OpenMP 回归比较完整稀疏结构、积分标签、CSF 对、贡献数和数值；
保留重复贡献，不先求和掩盖贡献丢失。现有 [mcp_compare](../implemented/mcp_compare.md)
可作语义比较工具，默认零容差仍不等于文件逐字节一致，二者分开记录。
比较工具目前以 `rcsf.inp` hash 检查身份；描述符输入实验应保留配套原文本基准，
以后才考虑增加经过验证的 canonical basis 身份协议，不用随意跳过身份检查代替证明。
同时让原下游 reader 读入，并在固定径向轨道下比较矩阵/RCI，不能只凭最终能量近似相等验收。

若浮点运算次序保持不变，首先要求精确一致。确需改变次序时，在运行前另定误差预算及理由，
不在看到差异后放宽阈值。所有通过结论只覆盖相应软件/数值行为，不建立活性空间收敛或物理精度。

## 7. 最初静态核查的边界

已做：读取上述当前源码、F2018/数据流/既有输入讨论；核对 OpenMP、GNU Fortran、
Parquet、Open MPI 和 Linux 官方规则；检查本文引用与 Git 改动边界。
没有完成深层角动量调用闭包的全量可重入审计，未复跑既有 MCP 工具测试，未实现或计时任何候选方案。

额外静态线索：`LODCSH2` 359–380 行的重复检查中，首次不相同会 `GO TO 17` 跳出整个检查区。
不能仅凭外层 `DO J=1,NCF-1` 就宣称它实际完整查重或据此认定加载存在二次复杂度瓶颈。
后续应以独立重复-CSF 夹具检查其契约；如确认缺陷，修复和预期差异单独记录，
不要混入新输入适配的等价性结论。本次未修改该逻辑，也未将其认定为现场性能根因。

原建议是先执行 B 的小型适配器与数组对照夹具，并建立 A 的阶段计时。
这部分现已实现；线程工作区、调度及单机样例已有首版实现和实测。显式 context、排序并行和更广性能矩阵仍待后续工作，详见文首说明。
