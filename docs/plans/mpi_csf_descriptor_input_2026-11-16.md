# MPI 描述符 CSF 输入接入与性能验证实施计划

日期：2026-11-16。源码核查基准：`caaf3a51ce12ec8e85e30622fc854c9e2a070683`。
状态：**代码分析及实施计划；本文新增功能尚未实施，也未运行新的构建、MPI/SCF 计算或性能实验。**
本轮开始时 suite 工作树干净。所有拟修改代码、测试和文档均归属 suite，不修改上游 `grasp/` 或 Python 子项目。

## 1. 结论与本次范围

当前源码与“MPI 都尚未接入描述符”的假设不完全一致：

| 程序/能力 | 当前源码状态 | 本计划的处理 |
| --- | --- | --- |
| `rangular` / `rangular_omp` | 已支持描述符二进制，保留文本入口。 | 作为共享 reader 的回归对象。 |
| `rangular_mpi` | **已支持** `--csf-descriptors`；rank 0 解码，广播头部和当前 block 数组；已有分块广播、计时和测试。 | 核对实际构建/调用的是此版本，复验既有功能；不再重复实现。 |
| `rmcdhf_mpi` | `SETMCP → SETCSLmpi → LODCSLmpiGG → LODCSH2GG`；CSF 入口仍固定为文本。 | 新增描述符入口。 |
| `rmcdhf_mem_mpi` | 相同的文本 CSF 加载链；MCP 在 `SETMCP` 中还会建立 memory 缓存。 | 与普通 MPI 版共用新增输入适配层；分别验证。 |
| RMCDHF 运行监测 | 已有外部 Linux 采样器，不等于 Fortran 输入阶段计时或新 reader。 | 保留采样器；本次只补充可关闭的输入阶段计时。 |

`rangular_mpi` 接入位于提交 `212bcc6`，且代码已在上述 HEAD 中：
[输入模块](../../src/appl/rangular_mpi/rangular_mpi_input.f90) 的
`parse_mpi_angular_options`、`load_mpi_angular_header`、`load_mpi_angular_block`，以及
[主程序](../../src/appl/rangular_mpi/genmcpmpi.f90) 87、139、168 行构成实际调用链。
[CMake](../../src/appl/rangular_mpi/CMakeLists.txt) 也已编译该模块。
源码存在不能证明服务器上安装的二进制或作业脚本已经启用该入口；本轮未核对服务器部署。

目标是形成以下可用链路：

```text
rCSFs Parquet + 匹配的 header.toml
  → suite prepare_csf_descriptors.py → 同一份 GCSDESC1 文件
  → rangular_mpi --csf-descriptors → 原 MPI MCP 分片
  → rmcdhf_mpi / rmcdhf_mem_mpi --csf-descriptors → 原 SCF 求解流程
```

两个 RMCDHF MPI 程序默认仍读取文本；新路径不再需要 `rcsf.inp`，但仍需要匹配的
MCP、`isodata`、`rwfn.inp` 及原有交互输入。这里的“读取描述符”指读取现有 suite 二进制交换文件，
**不是让 Fortran 直接解析 Parquet**。

本次不改角系数算法、MPI 列分工、MCP 格式、SCF/Davidson 算法、轨道改善顺序或收敛判据。
不新增求解器变体，不恢复已退役程序；串行 RMCDHF、RCI、CSFG 的新输入入口不属于本次交付。
普通 MPI 已移除的核极化入口保持现状，memory MPI 的可选核极化行为也保持现状。

## 2. 源码依据与必须保留的差异

以下行号对应文首基准；后续以链接和过程名为准。

### 2.1 rangular MPI 与 OpenMP 的差距不能归因于缺少描述符

- [rangular_mpi_input.f90](../../src/appl/rangular_mpi/rangular_mpi_input.f90) 102–175 行：
  rank 0 读取并装配，随后广播；描述符数组消息按 CSF 列分块，限制最大 `3*NNNW*列数`。
- [rangular_parallel.f90](../../src/appl/rangular/openmp/rangular_parallel.f90) 的
  `prepare_batch` 在当前 JA 列内分配 JB 任务，每批等待完成；`load_parallel_pair` 将结果交回主线程。
- [串行/OpenMP MCP 主流程](../../src/appl/rangular/mcp_gg.f90) 继续串行过滤、编号、累计和排序；
  [MPI MCP 主流程](../../src/appl/rangular_mpi/mcpmpi_gg.f90) 则由各 rank 分别处理局部列、排序并写分片。

这支持“当前实现的并行覆盖范围不同”，不支持“OpenMP 技术必然慢于 MPI”。
已有 [x86 诊断报告](../reviews/rangular_mpi_openmp_x86_diagnosis_2026-11-11.md)
记录 23,459-CSF 样本的 8 线程/8 ranks 中位耗时约 27.06/13.49 秒；
后续 [MPI 描述符报告](../reviews/rangular_mpi_descriptor_benchmark_2026-11-11.md)
记录同一 MPI 程序文本/描述符输入及广播约 0.411/0.061 秒，整程约 12.565/12.162 秒。
**这些是仓库保存的历史测量，本轮未复测；不同批次的整程时间不能直接拼接计算加速比。**

因此优先补齐 MPI 求解链；OpenMP 跨列调度、结果处理和排序优化应另立实验，不能把它们混入输入改造。

### 2.2 RMCDHF 不是 rangular 的逐 block 三数组加载模型

| 入口/状态 | 当前行为与实施约束 |
| --- | --- |
| [rmcdhf_mpi/rscfmpivu.f90](../../src/appl/rmcdhf_mpi/rscfmpivu.f90) 140–143 行；[memory 主程序](../../src/appl/rmcdhf_mem_mpi/rscfmpivu.f90) 142–145 行 | 先调用 `SETMCP`，再读 CSF，最后比较 NCORE。必须保留 MCP 先行及两类求解器各自的后续初始化。 |
| [SETMCP](../../src/appl/rmcdhf_mpi/setmcp.f90) 78–94、121–137 行 | 从 MCP 得到 NCORE、NBLOCK、NCFBLK、IDBLK、NELEC、NCF、NW、DIAG/ICCUT/LFORDR；检查 rank/进程数。`NCFBLK` 分配为 **0:NBLOCK**，且第 0 项为 0。 |
| [memory SETMCP](../../src/appl/rmcdhf_mem_mpi/setmcp.f90) | 还调用 `SETMCP_MEM`；这部分成本属于 MCP 输入，不能算作 CSF 解析。 |
| [SETCSLmpi](../../src/appl/rmcdhf_mpi/setcslmpi.f90) 57–95 行 | rank 0 打开文本、调用 `LODCSH`，广播轨道表；各 rank 分配全基组 `IQA(NNNW,NCFTOT)`。 |
| [LODCSLmpiGG](../../src/appl/rmcdhf_mpi/lodcslmpiGG.f90) 45–57 行 | rank 0 一次读全部 blocks，再广播 IQA 和 JPGG；现有 `NNNW*NCF` 单次 count 有默认整数乘法溢出的规模风险。 |
| [LODCSH2GG](../../src/appl/rmcdhf_mpi/lodcsh2GG.f90) 135–138、364–367、428 行附近 | 按文本原 block/CSF 顺序拼接全部占据数；每块首个 CSF 设置 `JPGG=(2J+1)*parity`，最后检查总 CSF 数。 |
| [公共描述符 reader](../../src/lib/lib9290/csf_descriptor_input.f90) | 已有格式、物理态和耦合校验；`install_descriptor_header` 会设置 NCF/NBLOCK，并重新分配 `NCFBLK` 为 1 起点布局；`read_descriptor_block` 要求三数组形状与该 block 完全一致。 |

**不能在 `SETMCP` 后直接调用 `install_descriptor_header`**：这会覆盖 MCP 给出的契约、
重新分配已有指针并丢失 `NCFBLK(0)` 语义。也不能临时把全局 `NCF` 改成当前块长度后忘记恢复。
RMCDHF 的全局 `NCF` 在新 loader 中始终表示全基组 CSF 数。

MPI 两个 `lodcslmpiGG.f90` 当前相同；两个 `setcslmpi.f90` 的主要差异是旧错误路径使用
`STOP` 或 `ERROR STOP`。适配器应共享一次实现，避免再出现两个独立描述符解码器。

### 2.3 保留低内存状态，同时处理遗留 EAL 依赖

普通 RMCDHF 加载器已停止长期分配 `JQSA/JCUPA`；主路径使用全基组 IQA 和每块 JPGG。
在 8 位 byte 条件下，CSF 主数组约为每 rank `NNNW*NCFTOT` 字节，
不应接入后变成 rangular 的约 `5*NNNW*NCFTOT` 字节。
这也不表示消除了 rank 级副本、MCP 缓存或 SCF 工作区。

仍有一个静态风险：[GETALDWT](../../src/appl/rmcdhf_mpi/getaldwt.f90) 的 `CASE(5)`
调用 [ITJPO](../../src/appl/rmcdhf_mpi/itjpo.f90)，后者读取未由当前 loader 分配的 JCUPA；
memory MPI 也保留此链。普通 rangular 输出的 `DIAG=.FALSE.` 通常进入 EOL，不能据常用路径成功就忽略 EAL。
本轮未运行触发实验，因此只确认依赖冲突，不宣称已复现崩溃。

实施时单独建立权重测试，再将这两个程序的标准权重来源改为按 block 展开的 `abs(JPGG)`，
保留权重归一化及其它选项。使用 `NCFBLK` 的局部前缀和定位全局 CSF，
**不能依赖此时尚未在 SCF 中分配的 `NCFPAST`**。
此项作为独立正确性修复记录，不以描述符加速掩盖；不修改通用 ITJPO 或恢复全量 JCUPA。
这项修复及输入等价性不构成完整 EAL 求解流程的数值验收。

## 3. 选定的实施设计

### 3.1 公共 reader：复用语义校验，新增仅输出占据数的入口

修改 [csf_descriptor_input.f90](../../src/lib/lib9290/csf_descriptor_input.f90)：

1. 保留 GCSDESC1 v1、原三数组接口和普通态支持范围，不新增格式版本或独立 seniority 表。
2. 新增 `read_descriptor_occupations_block(source, block_index, occupations, status, message)`。
   `occupations` 直接指向全基组 IQA 的连续 block 列区段；每块只读一次，维持 `next_block_index` 检查。
3. 将已有批量读取/逐行解码整理为共享内部实现。新入口仍经过同一个 `decode_descriptor_row`，
   验证占据、seniority、子壳层态、耦合、总 J、宇称和电子数；不能只取每四字节中的占据数。
4. 耦合及子壳层量子数只使用单行临时数组，解码后丢弃；原始描述符批缓冲沿用至多 8192 行。
   不分配全量或整个大 block 的 JQSA/JCUPA，不先装配三数组再复制 IQA。
5. `open_descriptor_source` 及新入口提供可返回的错误状态；所有可达读取、分配、形状、物理态错误均传播到调用者。
   为现有串行/rangular 调用保留原失败即终止的兼容包装，避免静默忽略错误。
   `lib9290` 内不引入 MPI 依赖。
6. 失败后关闭流并回收临时资源；测试同一进程内连续打开、失败、重新打开。

目前 reader 的 `require_descriptor`、`require_state`、`require_row` 会直接 `error stop`。
只在最外层增加一个 status 参数不能解决深层终止；状态传播必须覆盖整个调用链。

### 3.2 共享 RMCDHF MPI 适配器与应用入口

新增 `src/appl/rmcdhf_common/`，以一个小型静态支持目标 `rmcdhf_mpi_input_support`
提供 `rmcdhf_mpi_input.f90`，包含命令行选项、描述符/MCP 对照、MPI 加载及输入阶段计时。
两个 MPI 应用链接同一支持目标。它依赖现有 `mod/9290/mpiu90`，不复制数值库；
在根 CMake 的 `if(TARGET mpiu90)` 中、两个应用目录之前创建。
串行-only 配置不构建该目标，也不需要 MPI。

两个 `rscfmpivu.f90` 的改动相同但保持各自后续流程：

1. 在打开输出和启动 MPI 前解析参数，支持无副作用的 `--help`；拒绝未知选项、缺值和重复输入选择。
2. `startmpi2` 后按原规则确定 rank 工作目录与 rank 0 的 serial I/O 目录 `permdir`。
3. 保留原 `SETMCP`，保存其给出的完整元数据快照。
4. 描述符模式调用共享 loader；文本模式保留原 `SETCSLmpi` 和交互顺序。
5. 新 loader 成功后进入原 `GETSCDmpi → SCFmpi`，在 `stopmpi2` 之前关闭计时文件。

计划 CLI：

| 参数 | 契约 |
| --- | --- |
| 不传输入参数 | 原 `permdir/rcsf.inp` 路径，默认行为不变。 |
| `--csf-text FILE` | 显式文本路径。 |
| `--csf-descriptors FILE` | GCSDESC1 文件；与文本参数互斥。失败时不得偷偷回退文本。 |
| `--timings NEW_FILE.csv` | 默认关闭；每 rank 写 `NEW_FILE.rankNNN.csv`，拒绝覆盖已有文件。 |

相对 CSF 路径以 `permdir` 为基准，不能以执行 `chdir` 后的 rank 目录解释。
新输入路径用可分配字符串，不再经过主程序现有的 128 字符 `file_rcsl`；
旧 `CPATH/disks/MPI_TMP` 的工作目录长度限制仍存在，本文不扩大该接口。
输入/计时路径要防止与 MCP、临时文件、日志、sum、mix、rwfn 等输入输出冲突；
计时路径不能与选定 CSF、`isodata`、`rwfn.inp` 重合。
沿用现有 rank 文件命名和有效进程规模，不额外扩展 rank 编号格式。

### 3.3 描述符加载的顺序与状态所有权

```text
所有 rank: SETMCP → 保存本 rank 的 MCP 元数据
rank 0: 打开 GCSDESC1 → 验证格式及头部 → 与 MCP 元数据逐项比较
所有 rank: 交换成功/失败状态 → 对照本 rank MCP 与 root 元数据
所有 rank: 保留 NCFBLK(0:NBLOCK)、IDBLK、NCF → 安装/广播轨道信息
所有 rank: 分配 IQA(NNNW,NCFTOT)，同步分配结果
按原 block 顺序:
    rank 0: 验证/解码到 IQA(:,first:last)
    所有 rank: 同步状态 → 按完整 CSF 列分块广播 IQA 区段
所有 rank: JPGG(block)=(block_two_j+1)*block_parity
rank 0: 关闭输入；所有 rank: 保持原 SCF 所需状态 → GETSCDmpi
```

对照字段包括 NCORE、NW、NELEC、总 NCF、NBLOCK、每块长度及 IDBLK 的 J/宇称，
并验证 `sum(NCFBLK(1:NBLOCK)) == NCF`、第 0 项为 0。
用现有 `descriptor_block_label` 生成同格式标签比较；不重排、合并同对称性的 blocks，也不改变块内行序。
保存并保留 MCP 给出的 DIAG/ICCUT/LFORDR，不根据新文件重新推导相互作用范围。

轨道数组 `NP/NAK/NKL/NKJ/NH` 和电子数来自通过验证的 descriptor header；
不要调用会重建 block 状态的 `install_descriptor_header`。
JPGG 的有效范围为 `1:NBLOCK`，其余项初始化为 0，避免广播未定义尾部。
rank 0 为唯一打开二进制的进程；工作 rank 不要求能访问输入文件所在文件系统。

### 3.4 MPI 大数组与错误协议

在 `src/lib/mpi90/` 新增小型 CSF 通信辅助模块，提供 byte 占据数组的分块广播及输入阶段错误同步。
参考既有 rangular 实现和 suite 的状态同步模式，不复用带核极化语义的过程。
首版供 RMCDHF 使用；现有 rangular MPI 入口保留，避免把已经可用的入口重写作为接入前提。

- 尺寸乘法、总行数及 block 前缀和用 `int64` 计算；进入 default-integer 数组/内核前检查范围。
- 消息按完整 CSF 列切片，`NNNW*columns <= huge(default integer)` 后才转换 MPI count；
  使用与 `integer(BYTE)` 宽度匹配的类型，GNU 15 目标上核验 `MPI_INTEGER1` 与 8 位存储。
- 每个区段都在最终 IQA 上连续，避免隐式生成全数组临时副本；在 Check 构建检查切片及临时数组诊断。
- 测试接口可传很小的 count 上限，强制多次广播；用宽整数算术测试极大逻辑尺寸，无需分配数 GB 夹具。
- root 读取/解码状态在 payload collective 前广播；所有 rank 的分配/本地契约错误先归约定位失败 rank。
  出错输出阶段、block/行和原因，然后 `MPI_Abort` 非零退出；不能让其它 rank 永久等在广播。
- MPI 调用自身失败直接按统一错误路径中止，不在失效 communicator 上继续复杂恢复。

新增状态 API 只消除新输入链对单 rank `ERROR STOP` 的依赖；旧 `SETMCP`、文本 parser 等仍有历史
`STOP/ERROR STOP` 分支。逐项测试本次涉及的失败路径；不能宣称全部 suite 错误退出已统一。
将两个 `lodcslmpiGG.f90` 的文本 IQA 广播也切换为同一分块 helper，保持文本/描述符具有同样的 count 保护。
这不等于扩大 MCP 位置号、矩阵索引或其它旧内核整数容量。

### 3.5 MCP 身份边界

现有 MCP 头没有完整轨道表及逐 CSF 身份摘要，GCSDESC1 的 `header_sha256` 也只绑定 header。
所以“计数和 J/宇称相同”不能证明 MCP 属于该描述符，更不能检测同一块内的 CSF 置换。

本期保留 MCP 格式，明确两层检查：

1. Fortran loader 强制检查所有现有可比元数据及 rank 配置，阻止可检测的不一致。
2. 验证/性能运行器在隔离目录中用同一份固定输入生成 MCP，再复制原始分片给各次 SCF，
   保存 descriptor、可选对照文本、MCP 全部文件、可执行文件、交互输入的 SHA-256 和来源信息。
   这只证明该受控运行链的对应关系，不把任意历史 MCP 重新命名后视为配对成功。

生产输入若缺乏可核验的 MCP 来源，应从选定基组重新生成对应 rank 数的 MCP。
新 MCP 身份 sidecar/格式升级可作为后续独立协议，本次不虚构现有 loader 能检查不存在的身份字段。
`rangular_omp` 的单套 MCP 不能改名为 MPI 分片后使用；进程数改变时必须重新生成正确分片。

## 4. 文件级工作清单

下表中的新路径、过程及测试名均为**计划新增**，不是当前已存在接口。

| 文件/位置 | 具体工作 |
| --- | --- |
| `src/lib/lib9290/csf_descriptor_input.f90` | 共享批量解码；占据数 block API；可返回状态的内部错误链及兼容包装。 |
| 新 `src/lib/mpi90/csf_input_mpi.f90`、其 CMake source list | 有界 IQA 广播、状态同步及失败中止；不引入第二份库。 |
| 新 `src/appl/rmcdhf_common/rmcdhf_mpi_input.f90`、`CMakeLists.txt` | 两程序共用选项、MCP 对照、轨道安装、全基组加载和最小输入计时。 |
| 根 `CMakeLists.txt`；两个 MPI 应用的 `CMakeLists.txt` | MPI 条件内增加支持目标及依赖；使用现有模块目录机制，避免并行构建 `.mod` 冲突。 |
| 两个 `rscfmpivu.f90` | 参数解析、选择加载入口、保留 MCP 快照、输入计时；保留各自 GETSCD 接口。 |
| 两个 `lodcslmpiGG.f90` | 原文本解析保留；替换可能溢出的单次 IQA 广播。 |
| 两个 `getaldwt.f90` | 独立修复标准权重的 JCUPA 依赖；局部前缀和 + JPGG，并增加权重单测。 |
| 新 `test/rmcdhf/csf_input_arrays.f90`、`test_mpi_descriptors.py`；`test/CMakeLists.txt` | 数组和元数据对照、失败退出、输入路径、真实 MPI 链路和两求解器回归。 |
| 新 `scripts/compare_rmcdhf_mpi_inputs.py` | 同一求解器文本/描述符交替重复运行，保存输入/MCP/程序身份、计时、输出及明确的成功判定。 |
| 现有 `test/test_csf_descriptors.py`、`test/rangular/test_mpi_descriptors.py` | 覆盖公共 reader 变更后的旧接口；不削弱已有精确比较。 |
| 文档索引、实施说明、脚本 README、F2018 清单/警告基线 | 完成后登记支持范围、命令和实测证据；按实际新编译单元更新清单，不盲目提高警告预算。 |

若维护旧 Makefile/BUILDCONF 入口，其依赖清单要同步；支持的构建验收仍使用 CMake。
不在两个应用目录各复制一份共享支持模块。

## 5. 实施顺序与完成条件

| 阶段 | 工作 | 进入下一阶段的条件 |
| --- | --- | --- |
| P0 基线 | 固定 HEAD、工具链、库、构建选项和原始程序；确认实际 rangular MPI 支持选项；建立同 rank 数 MCP 夹具及两个 RMCDHF 的文本基线。 | 基线可重复；成功标记、实际迭代数、停止原因和文件完整性均记录，未把退出码 0 单独当成功。 |
| P1 reader | 新增占据数入口和可返回错误；用旧三数组 reader/文本 GG reader 核对有效数组。独立复现并修复标准权重依赖。 | IQA、JPGG、元数据精确一致；不分配全基组耦合数组；错误、重开及权重测试通过。 |
| P2 MPI | 实现 count 安全广播和状态协议；文本 GG 路径改用同一 helper。 | 1/2/4 ranks、强制小 chunk、非 root 分配错误、root 数据错误均无挂起；文本结果保持一致。 |
| P3 两应用接入 | 加 CLI、MCP 对照和轻量输入计时；同时接入普通/memory MPI。 | 描述符目录没有 `rcsf.inp` 仍能完成规定 SCF 步数；旧默认入口和两程序各自语义保留。 |
| P4 数值/构建回归 | 同输入类型、同 ranks 对照；公共 reader 的 rangular/OMP/MPI 回归；ordinary Release/Check 全量检查。 | 内部状态和科学输出符合第 6 节要求；无新增运行时错误，受影响配置可编译；已知限制单独记录。 |
| P5 性能与交付 | 执行受控输入/整程计时，更新实施说明、对照报告与运行示例。 | 分别报告加载与整程收益、预处理成本和内存；性能改善与否都如实记录，不能以预期代替结果。 |

每阶段保持可独立检查的变更；只有在用户要求时提交 Git。完成 P4 可确认软件输入能力，
P5 才能对指定体系和机器作性能结论，二者都不建立活性空间收敛或物理精度。

## 6. 验证矩阵与验收标准

### 6.1 loader 与状态测试

- 对两个实际文本 `SETCSLmpi/LODCSH2GG` 建立驱动，对比 descriptor 的
  NCORE、NELEC、NW、轨道顺序、总 NCF、`NCFBLK(0:NBLOCK)`、IDBLK、有效 JPGG 和全部 IQA。
  额外检查输入之后 block 指针下界及 MCP 元数据没有被覆盖。
- 覆盖无 core/有 core、空满壳、单/多开放壳、部分占据 J=0、显式 seniority、正负宇称、整数/半整数 J、
  多 block、同 J/宇称的独立 blocks，以及原始 CSF 行序。歧义/超支持态按现有契约拒绝。
- 对标准权重使用解析期望：一个 block 内所有 CSF 的未归一化权重为 `abs(JPGG(block))`；
  多 block、负宇称、块边界和归一化分别检查，不调用完整 EAL 求解器充当唯一 oracle。
- 非法 magic/version/端序标记、错误长度、缺失文件、截断、非法态/耦合、容量超限、block 元数据不匹配都应明确失败。
  匹配计数但行置换的输入不应宣称可由 MCP 头检测；在受控运行器层验证其来源被区分。
- 原三数组接口和新占据数接口共用物理语义，禁止以简化校验换取通过率或速度。

### 6.2 MPI 集成与真实小型 SCF 回归

计划新增 CTest 名：`rmcdhf_csf_input_arrays`、`rmcdhf_mpi_descriptors`、
`rmcdhf_mem_mpi_descriptors`、`rmcdhf_csf_mpi_errors`；现有 rangular 测试名保留。

测试矩阵至少覆盖：

| 维度 | 要求 |
| --- | --- |
| 求解器 | 两个 MPI 版各自与自身文本基线比；不得用普通版与 memory 版互为唯一数值基线。 |
| rank 数 | 1/2/4，包含没有分配到计算列的 rank；每种 rank 数重新生成匹配 MCP。 |
| 输入 | 默认文本、显式文本、描述符；描述符用例在目录中完全不放文本文件。 |
| 交互/模型 | full interaction 和一阶 ICCUT；小型多 block、多 CSF EOL、冻结轨道及非冻结轨道；固定一次与多次 SCF 更新。 |
| 路径 | 绝对/相对、含空格、长输入文件路径、rank 0 独占读取、timings 防覆盖与参数错误。工作目录遵循旧 CPATH 限制。 |
| 错误 | root header/payload 错误、非 root 分配/元数据错误、不同 rank 配置 MCP；检查非零退出和超时/残留进程。 |
| memory/NP | memory 版缓存路径确实被执行；关闭 NP 的主基线，以及 memory 版可选 NP 路径不受影响的现有专项回归。 |

先保存同构建、同 rank 数旧文本/新文本/新描述符结果，分别识别共享改动及输入表示的影响。
要求 IQA/JPGG 和所有离散元数据精确一致；输入改造不改变浮点操作顺序，
因此同程序、同配置首先要求 `rwfn.out`、`rmix.out` 及科学数值记录精确一致。
日志只剔除明确列举的时间、路径和输入格式诊断，不能过滤警告、能量、迭代数或停止消息。
多次 MPI 运行若基线本身有波动，先记录并定位，不能看到差异后再放宽阈值。
MCP 分片在复制入运行目录前记录完整文件集和 hash，避免 memory/SCF 运行改变文件状态后污染下一轮。

现有 [compare_first_sweep.py](../../test/f2018/compare_first_sweep.py)
CLI 仅支持串行 `rmcdhf/rmcdhf_mem`；可复用其输出完整性、非有限数和停止状态检查思想，
**不能直接把它当作现成的 MPI 端到端验收器**。
新的 MPI 驱动必须自行管理 launcher、rank 目录、MCP 分片、超时和输入选项。

### 6.3 构建与工具链

以下为**实施后**从 suite 根目录执行的 fish 示例，本轮未执行：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
/usr/bin/gfortran-15 -dumpversion
mpifort --showme:command
set suite_python (path normalize "$PWD/../graspkit-tools/.venv/bin/python")

cmake -S . -B build-csf-input-release \
    -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
    -DCMAKE_BUILD_TYPE=Release -DGRASP_FORTRAN_PROFILE=f2018-release \
    -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
    -DGRASP_BUILD_CSFG=OFF -DGRASP_RANGULAR_OPENMP=ON \
    -DPython3_EXECUTABLE="$suite_python"
cmake --build build-csf-input-release -j4
ctest --test-dir build-csf-input-release -N
ctest --test-dir build-csf-input-release \
    -R 'csf_descriptor|rangular_mpi_descriptors|rangular_openmp|rmcdhf.*(descriptor|csf)' \
    --output-on-failure

cmake -S . -B build-csf-input-check \
    -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
    -DCMAKE_BUILD_TYPE=Debug -DGRASP_FORTRAN_PROFILE=f2018-check \
    -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
    -DGRASP_BUILD_CSFG=OFF -DGRASP_RANGULAR_OPENMP=ON \
    -DPython3_EXECUTABLE="$suite_python"
cmake --build build-csf-input-check -j4
ctest --test-dir build-csf-input-check --output-on-failure
ctest --test-dir build-csf-input-release --output-on-failure
```

先检查 `ctest -N` 确实包含新增测试，不能让过滤表达式匹配零项后报告通过。
补充 MPI-only、serial-only、关闭 OpenMP 的受影响构建配置；公共 9290 改动还要覆盖已有属性程序测试。
更新 F2018 source/target 清单并执行相应维护检查；记录缺失依赖、平台限制和基线失败。
不以跳过用例或降低检查等级来获得表面全绿。

## 7. 性能实验设计

**目的与假设：** 在科学流程完全相同的条件下，检验二进制描述符能否减少 RMCDHF 初始化的
文本解析/装配时间；预期输入阶段可能改善，整程收益取决于其占比，数值和幅度目前未知。
此改造不预期直接解决 SCF 全程反复 D 状态、MCP 复读、Davidson 通信或轨道改善成本。

新增输入计时沿用已有 CSV 列约定，在 RMCDHF 支持模块内实现最小功能；
不让 RMCDHF 链接 rangular 应用目标，也不为本项迁移现有 rangular profile 模块。

| 阶段 | 边界/含义 |
| --- | --- |
| `mcp_setup` | 原 SETMCP，memory 版包含 SETMCP_MEM；单独报告。 |
| `csf_input_total` | 整个文本/描述符适配；用于公平的主输入对照。 |
| descriptor 子阶段 | root 头部/元数据检查、各 rank 分配、root 解码、payload 广播分别记录；子项嵌套于 input total，不重复相加。 |
| `program_before_finalize_including_interaction_input` | startmpi2 返回后至 stopmpi2 前；必须使用固定 stdin，并另报进程外完整墙钟。 |

不为计时新增 barrier，不逐 CSF 记录日志。通信时间包含等待，不能称为纯网络传输；
各 rank 时间彼此重叠，不能把阶段最大值之和当成整程耗时。

实验步骤：

1. 固定程序 hash、GNU 15/MPI/BLAS、绑定/核数、数学库单线程、初始轨道、核/网格、EOL 根/权重、
   冻结轨道、交互 stdin 和迭代控制；普通与 memory MPI 分开测量。
2. 先在小型夹具完成正确性验收，再选实际有代表性的 CSF 规模。每种 ranks 使用匹配的原始 MCP 副本。
   对同一程序交替文本/描述符，至少三轮，保存每次原始耗时及中位数/范围。
3. 输入复制、Parquet 预处理、输出校验及归档在计算计时之外；一次性转换成本另外报告。
   既报告已有二进制复用，也报告转换加运行的端到端成本。记录缓存/文件系统条件，不把默认临时目录当生产磁盘。
4. 同时报各 rank 峰值 RSS、总作业内存、读写量；观察 root 的临时内存，验证没有恢复全量 JQSA/JCUPA。
   外部采样按现有工具执行，注明采样粒度不能保证捕获瞬时峰值；检查开关计时/采样的扰动。
5. 计算 `input_fraction=T_input_text/T_total_text`、输入节省及整程节省。
   在其它工作不变的假设下，最多只能消除原输入耗时；如整程差异明显超出这一量级，先检查噪声、迭代轨迹和 I/O 条件。

不预先设定未经证据支持的加速倍数。若加载变快但整程几乎不变，软件接入可通过，
性能结论应明确写出收益范围；若整程变慢则分析原因，不将新模式强制设为默认。
实际科学实验开始前需在会话说明具体输入、目标、对照、预期和判据；本计划没有启动实验，也不授权 SLURM 提交。

## 8. 本轮完成与后续交付

本轮完成了入口、共享 loader、MPI 广播、MCP 契约、CSF 状态消费者、构建及测试设施的静态核查，
并核对已存在的 rangular MPI 描述符实现与历史性能报告。
只新增本计划并更新文档索引；没有修改生产代码、重建程序、重跑既有测试或产生新的性能/物理结果。

实施完成后应交付：两程序的可选描述符入口、共享占据数 loader、受控失败及 count 安全广播、
数组/MPI/SCF 回归、可重复对照运行器、当前使用说明，以及明确区分软件验证和性能测量的结果报告。
下一步从 P0 开始，先确认实际使用的 rangular MPI 构建身份并冻结两个 RMCDHF 文本基线。
