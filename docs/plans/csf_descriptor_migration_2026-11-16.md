# 任务二：全部普通程序迁移至 Rust 原生 CSF 描述符输入

创建：2026-11-16；方案修订：2026-11-17。状态：**生产输入迁移与功能回归已完成；既有历史告警门禁失配单列，不宣称其通过。**
当前输入与构建契约见[实施说明](../implemented/native_csf_descriptors.md)，
验证证据和既有告警门禁边界见[验收记录](../reviews/native_csf_descriptor_validation_2026-11-17.md)。

本轮源码基准：suite `5806884a6c448fab7c6bc8ca4ccfe84908bbd5b8`；
rCSFs `fea061fc741c40fe47ad00ead36ad20d55d34398`（Cargo 包版本 3.0.1）。
最初交接基准为 `6f96222`，旧输入诊断证据基于 `9d6ee31`，不能用旧状态代替实施时复核。

用户已明确最终要求：**suite 中全部保留的 CSF 消费程序只读取描述符；由进程内 Rust 库直接读取
Parquet 和配套 header TOML，经 C ABI 交给 Fortran，不依赖 Python 转换、GCSDESC1 中间文件或旧文本输入入口。**
这取代原计划的 `prepare_csf_descriptors.py → GCSDESC1` 生产数据路径。
本计划已实施；改动按 rCSFs、suite 和父仓库 gitlink 分层提交，未推送。另见[重命名任务](variable_renaming_followup_2026-11-16.md)和
[静态诊断任务](static_diagnostics_fixes_2026-11-16.md)。

## 1. 实施前现状与需要阅读的资料

当前保留的 `rangular`、`rangular_mpi` 已支持 **GCSDESC1**，尚不能直接读取 Parquet；`rwfnestimate` 已支持 header TOML。
四个普通 RMCDHF、RCI 及其他性质程序仍需逐入口盘点，不能将 rangular 支持误认为全 suite 支持。
整理期间另一任务已将 rangular OpenMP 移除提交为 `6f96222`；旧文档里的 `rangular_omp` 仅作历史参考，不恢复该目标。

必读资料：

- [现有描述符契约](../implemented/rangular_descriptor_input.md)、[MPI 输入契约](../implemented/rangular_mpi_descriptor_input.md)。
- [早期 RMCDHF MPI 细化计划](mpi_csf_descriptor_input_2026-11-16.md)：低内存布局、MCP 对照、广播与验证很重要；
  其中旧变量名、OpenMP 目标、命令、默认入口策略、GCSDESC1 路径及实施状态须按当前源码和本文校正。
- [rwfnestimate header 输入](../implemented/rwfnestimate_header_input.md)。
- [剩余非风格诊断](../reviews/remaining_nonstyle_static_analysis_2026-11-16.md)：CSF 错误输入的实际复现和合法闭核边界。

最终生产数据路径（待实施）：

```text
rCSFs V2 描述符 Parquet + 字节匹配的 header TOML
  → rCSFs 所有的无 Python Rust reader（分批 Parquet 解码、TOML/格式/绑定校验）
  → 版本化 C ABI（同一进程，调用方缓冲区）
  → suite 公共 Fortran 适配层（物理校验、J/K 补全和数组装配）
  → 调用方需要的原计算数组与元数据
  → 现有数值内核
```

V2 文件契约、C ABI 版本和历史 GCSDESC1 版本分别管理，不能混淆。Rust 不生成临时 `.c` 或
GCSDESC1 文件，不启动 Python、rCSFs CLI 或其他转换子进程；C ABI 不要求编写额外的 C 程序。

初期及本任务验收使用 **Parquet + header TOML** 双文件契约。header 内的五行字符串是现有描述符配套元数据，
解析这些元数据不等于保留 `.c`/`rcsf.inp` 的三行 CSF reader。单文件 Parquet 是另一项协议变更，
需将 core、轨道顺序、显式 block lengths 等必需信息完整嵌入并版本化，不能假设当前 footer 已含全部数据。

### 1.1 实施前 rCSFs 可复用代码盘点

以下路径相对工作区的 `rCSFs/`，链接定位本轮实际读取的源码；“可复用”不等于已完成新接口验证。

| 现有代码 | 可复用内容 | 接入限制/必须补齐 |
| --- | --- | --- |
| [descriptor_schema.rs](../../../rCSFs/src/descriptor_schema.rs) | `DescriptorVersion`、`DescriptorLayout`、`MISSING`、角动量上界、`validate_descriptor_storage`、`validate_seniority_encoding`、`validate_output_schema`、`hash_header_file` | 作为 V2 契约唯一来源；补统一、严格的输入元数据校验，不能只复用其中两个标记检查就宣称完整兼容 |
| [csfs_descriptor.rs](../../../rCSFs/src/csfs_descriptor.rs) 的 `restore_v2_descriptor_parquet_to_outputs` | `ParquetRecordBatchReaderBuilder`、footer 读取、Int8/空值检查、批次遍历、显式 block lengths 及计数核对 | 提取 reader，不调用此恢复函数：它创建文本输出；其 hash 校验是“存在时才检查”，新生产入口必须要求存在；部分键使用首个匹配值，需拒绝重复/冲突键 |
| 同文件的 `parse_restore_header_from_toml` 与 `parse_peel_subshells_from_header_lines` | 现有 TOML 字段及 peel 顺序解析思路 | 前者为私有恢复辅助函数，后者为 `pub(crate)` 且拒绝空 peel、会过滤部分 token；新 header API 须严格校验 core/peel、标签、重复轨道并支持纯闭核，不能原样暴露 |
| [csfs_conversion.rs](../../../rCSFs/src/csfs_conversion.rs) 的 `HeaderInfo`、`BlockInfo`、`HeaderData` | 已有 serde 数据模型；电子数是 `Option<usize>`；显式 block lengths/count 和 conversion stats | 分离读取所需元数据与转换统计，核对存在的冗余计数；使用宽整数和 checked arithmetic；旧 header 缺 electron_count 的处理见 §2.2 |
| [descriptor_v2.rs](../../../rCSFs/src/descriptor_v2.rs) 的 `decode_v2_into` | 一行解码到复用缓冲，保留 `-1` 与显式零，物理 seniority/打印标签转换 | 可作格式回归和解码参考，但不是完整 suite 物理校验器；输出是恢复文本用的 occupied/state 结构，不直接作为 Fortran ABI，不绕道打印标签再重建物理 v |
| `descriptor_v2.rs::restore_file` 及 `csfs_descriptor.rs::restore_csfs_from_descriptors(indices=Some(...))` | 原顺序/显式块边界的测试参照 | 需要全量 rows 或 `CompleteCsfFile`，不作为低内存生产 reader；不能通过恢复文本再调用旧 Fortran reader |
| [seniority.rs](../../../rCSFs/src/seniority.rs)、[shared_generation](../../../rCSFs/src/configuration_generation/shared_generation/mod.rs) | 既有状态表、`Subshell` 标签及 n/kappa/capacity 解析 | 不复制状态表。Rust 与 suite 支持范围不完全相同；例如 Rust `Subshell` 约束 `n<=15`，状态表接受单电子 `2j<=21`，不能机械替代 suite 的范围检查 |
| [Cargo.toml](../../../rCSFs/Cargo.toml)、[lib.rs](../../../rCSFs/src/lib.rs)、[build.rs](../../../rCSFs/build.rs) | 已有 `rlib/staticlib/cdylib`；PyO3 为可选 feature，原生 CLI 有无 Python 构建路线 | 默认仍启用 `python`；算法模块大多无条件编译，关闭 Python 不等于 reader-only。当前 `src/` 未发现导出 C ABI；macOS build.rs 的 cdylib 参数也需按产物用途审查 |

格式复用与物理支持必须分开：`validate_record` 检查占据、局部态及耦合位置，不检查完整耦合三角链、
全电子数、总宇称等全部 suite 条件；它还拒绝空 occupied 集合。因此不能直接把它作为纯闭核或所有
suite 输入的唯一合法性判断。suite 仍使用共享 Fortran 物理适配器和 `terms_C`，Rust 的格式契约留在 rCSFs。
复用生成器类型时必须核对其范围限制；若会缩窄 suite 已支持的合法输入，应在 rCSFs 中分离通用读取类型与
生成器限制，而非静默缩窄范围。已有 Rust/Fortran 状态表的交集需做差分测试，不在本任务中扩张数值内核支持集。

可复用的测试来源包括 rCSFs 的 [descriptor_i8_test.rs](../../../rCSFs/tests/descriptor_i8_test.rs)、
[csfs_descriptor_test.rs](../../../rCSFs/tests/csfs_descriptor_test.rs)、
[native_cli_test.rs](../../../rCSFs/tests/native_cli_test.rs)、
[header_electron_count_test.rs](../../../rCSFs/tests/header_electron_count_test.rs)，
以及 `descriptor_schema`、`descriptor_v2`、`seniority` 的模块测试。
它们是已有测试定义，不是本轮新运行结果；rCSFs 维护测试仍须独立于 suite/外部 GRASP 二进制。

### 1.2 仓库归属与构建交付

- **rCSFs**：通用 Rust 描述符/header reader、严格 V2 文件契约、可选 C ABI、对应 Rust/FFI 测试。
  在所属仓库演进或抽取共享实现，不向 suite/父仓库复制一套 schema 或状态表。
- **suite**：`ISO_C_BINDING` 声明、内部数组装配、应用状态/MCP 对照、MPI、CMake 链接及集成测试。
- **Tools**：仅在流程/脚本需要调整时在其仓库改动；最终生成直接输入 Parquet/header 的调用，
  不能残留 Python 转 GCSDESC1 或恢复 `.c` 才能启动 suite 的步骤。Tools 本身的 Python 功能不因此删除。

优先在 rCSFs 增加可选 reader/FFI 构建边界，向 suite 提供 `staticlib` 和版本化 C 头文件，
明确关闭默认 `python`，不启用 CLI。模块和依赖按 reader 所需功能隔离；如果当前单 crate 无法形成
有效隔离，在 rCSFs 仓库内抽取共享核心/reader crate，现有生产者也依赖它，不复制实现。
具体新增 feature、crate、符号名称在 D1 定稿，本文不把它们写成已存在的构建命令。

构建需要 Rust/Cargo 和现有 GNU Fortran 15；Parquet/Arrow 使用 Rust 实现，不要求安装 Arrow C++。
锁定 Rust 工具链、Cargo.lock、启用的压缩功能和库来源 revision/version；CMake 不依赖开发者的可编辑 Python
安装，也不能隐式依赖相邻 checkout 的未提交状态。提供明确的本地开发源码覆盖和可发布依赖来源，
并验证依赖已备齐的离线构建。压缩后端可能仍需要 C 编译器；静态链接不等于所有平台均无系统动态依赖。
验收时检查最终链接依赖及无 Python 符号，在没有 Python/Rust 工具链的运行环境启动发布产物；
MPI、Fortran runtime 和系统库依赖照常记录。该构建/运行验证尚未执行。

### 1.3 应用覆盖表（基线盘点；当前均已接入，验收见实施说明）

| 保留的程序族 | 当前入口/特殊需求 | 最终验收 |
| --- | --- | --- |
| `rangular`、`rangular_mpi` | GCSDESC1 或文本；逐块三数组 | 改接 Rust，保持原块/CSF 顺序、ICCUT 及 MCP 输出 |
| `rmcdhf{,_mem,_mpi,_mem_mpi}` | 文本 `setcsl/setcslmpi`；SETMCP 先行、全基组占据 | 低内存加载、MCP 核对、标准权重依赖核查及必要修复；四个入口全部验收 |
| `rci`、`rci_mpi` | 文本头部扫描和 matrix 内逐块加载 | 头部与矩阵阶段共用同一描述符源及块顺序 |
| `rbiotransform`、`rbiotransform_mpi` | `setcslb/lodcslbio` 及应用内 setcsla 变体；两套基组 | 独立句柄、轨道映射、重读和变换后产物被下游直接消费 |
| `rtransition`、`rtransition_mpi`、`rtransition_phase` | `setcsl/setcslm`、`ldcsl1/ldcsl2` 及 MPI 变体；左右态 | 两套基组和全部读取/重读阶段均走描述符 |
| `gj{,_mpi}`、`rhfs{,_mpi}`、`ris4{,_mpi}`、`rdensity` | 主要经共享 `setcsla/lodcsl` | 各可执行目标和分支独立确认，不因共享函数迁移就跳过集成 |
| JJ2LSJ 的串行/OpenMP 构建 | `jj2lsj_code` 经 `setcsla`；标签/混合系数流程 | CSF 直接输入描述符，标签输出和 mixing 格式保持原契约 |
| `rwfnestimate` | 已有 header TOML 与文本入口 | 只保留描述符配套 header 输入，不扫描无关 CSF payload，删除旧 CSF 文本分支 |

此表依据当前 CMake 保留目标和调用入口建立，不代替 D0 的全调用图核查。继续检查共享库、辅助程序、
隐藏的再次 OPEN/REWIND、生成脚本和程序间产物；不消费 CSF 的程序不强加 Rust 读取依赖。
如果上游程序生成的 CSF 被下游 suite 再读，必须提供匹配 header 的描述符产物，不能保留内部文本绕行。

## 2. 公共接口与状态所有权

从现有 [csf_descriptor_input.f90](../../src/lib/lib9290/csf_descriptor_input.f90) 演进，
将 GCSDESC1 文件解码替换为 Rust C ABI 数据源，保留共享的物理语义校验和装配逻辑。
对外需要头部加载、逐块/分批读取、按需全块装载、关闭和受控错误返回。
新接口名及模块拆分是实施设计，本文不将拟新增路径描述成已经存在。

1. 分开“读取并校验数据”和“安装到应用全局状态”。现有 `install_descriptor_header` 会修改全局计数并分配块数组，
   不能无条件供所有程序调用；重复加载、失败后重开和上下态两套基组必须有明确所有权。
2. 需要三数组的应用继续获得占据、子壳层量子数、耦合数组；RMCDHF 增加只输出占据数的块入口。
   占据数入口仍验证全部描述符通道，耦合校验使用单行/有界批缓冲，不恢复全基组 JQSA/JCUPA。
3. 深层 `require_descriptor/require_state/require_row` 的失败要真正传播到调用方；只在最外层增加 status 参数不够。
   失败时关闭流、回收临时资源，不能让调用方使用部分安装的元数据或数组。
4. Rust reader 和 `lib9290` 保持无 MPI 依赖。MPI 模块负责 root 解码、元数据/数据广播、各 rank 分配状态归约和非零 `MPI_Abort`。
   先同步状态再进入 payload collective，按 CSF 列拆分 count；不依赖单 rank STOP 的启动器行为。
5. 保留原轨道、块、CSF 顺序；不合并相邻同 J/宇称块，不改变一阶参考区间。
   明确各应用的空轨道剔除、闭核、重排及左右态映射，不能照搬 rangular 的安装过程。
6. 支持的物理态、seniority、w、J/K 编码和整数容量沿用已验证契约；未支持的态明确拒绝，不静默补猜。
   64 位文件计数不等于旧内核支持任意大基组。

`rwfnestimate` 只消费轨道表与总电子数，可保留 header-only 能力；无需为统一输入而扫描全部 CSF 行。
纯闭核基组要有明确合法表示，不能靠失败的 READ 推测出一个 CSF。

### 2.1 Rust/C/Fortran 边界

- 提供打开、查询元数据/有序块表、分批读入调用方缓冲、关闭和错误消息接口，并有 ABI 版本检查。
  Rust 句柄持有文件、reader 和有限批缓冲；允许同时存在左右态两个句柄，不使用全局“当前输入”。
- Fortran 通过 `BIND(C)` 接收固定宽度整数和显式长度/容量；不跨边界传递 Rust `Vec`、`String`、
  enum、Arrow 对象或 Fortran 非互操作派生类型。计数使用 64 位，进入旧内核/default INTEGER 前检查范围。
- 首版 ABI 返回原始物理 seniority 的有序 Int8 通道；明确 `values(channel, row)` 的 Fortran 连续布局，
  Rust 零基偏移为 `channel + row_width * row`。不让打印标签恢复 API 成为计算数组的必经步骤。
- Arrow 列数据到行缓冲的打包允许有界复制，不宣称零拷贝。批次可以跨 row group，但对外块入口
  不越过显式 block 边界；不按 J/宇称聚类、排序或过滤行，禁止 materialize 整份 `Vec<Vec<i8>>`。
- 明确 EOF、有效行数、缓冲不足、格式错误、I/O 错误和资源错误；失败输出无效，不安装半成品状态。
  `close` 释放 Rust 资源，Fortran 释放自己的数组；初始化失败、部分读取失败、关闭后重开均可清理。
- FFI 不允许 Rust unwind 越过 C/Fortran；正常错误用 `Result` 转状态码，明确 panic 捕获和构建策略。
  不能声称 `catch_unwind` 能处理 abort/OOM。关键分配尽可能可失败返回，无法恢复的进程终止另列边界；
  MPI 可控分配错误必须进入状态归约并非零 Abort。

### 2.2 严格输入契约与元数据准备

Rust 新 reader 至少覆盖现有 Python 预处理的全部格式检查：从 **Parquet footer** 读取唯一且一致的
`descriptor_version=2`、`descriptor_dtype=int8`、`max_two_j=126`、`seniority_encoding=physical`、
`channels_per_subshell=4`、`subshell_count`、`missing_sentinel=-1`、`normalized=false`、
`peel_subshells`、`feature_columns`、`global_columns` 和必需的 `source_header_sha256`。
检查实际 Arrow schema 的列名/顺序/类型/nullability 和每批真实空值，不只相信 metadata；拒绝普通
三行 CSF Parquet、旧版本或归一化 ML 特征。不要把 Parquet 的物理 INT32 存储误判为逻辑 Int8 不合法。

header 字节只读一次并对同一快照计算 SHA-256 和解析 TOML，避免先解析再重开哈希产生配对竞态。
严格校验五行标识、core/peel 顺序与互斥、轨道范围、正 block lengths、宽整数求和、footer 总行数和存在的冗余计数。
`source_header_filename` 仅为提示，不替代内容绑定。所有文件句柄为只读，不重写 header。

现有 header 主要提供块长度，**不能假设已提供逐块 J/宇称和必需的 electron_count**。
需要完整块表的调用方在安装状态前进行有界预扫描：从原行序获取逐块对称性及全电子数，
检查每行占据/电子数和块内对称性一致，补齐元数据后重置同一输入源，再逐块完成所有物理通道验证。
不落地中间文件、不缓存全基组；重复读取成本纳入计时。若以后新增可靠的摘要元数据，须另立契约后才能省去预扫描。
文件在加载期间不得被修改，检测到不一致即失败；不把 header 哈希当作 payload 身份保证。

旧 header 缺电子数时，完整 reader 可从描述符行确定并校验；header-only 的 `rwfnestimate` 必须要求
明确、合法的含 core 电子数，缺失时报错，不能转去读取 `.c`。header-only 只证明头部合法，不证明 payload 合法。
纯闭核表示为非空 core、空 peel、一个长度为 1 的块及 `[total_two_j=0, parity=+1]` 描述符行；
空 peel 不能触发“缺数据时猜一个 CSF”。需在 rCSFs reader/相关生产者及 suite 夹具中共同验收。

## 3. RMCDHF 的关键约束

- 保留 `SETMCP` 先行及普通/memory 各自初始化流程，先保存 MCP 元数据再与描述符核对。
- `SETMCP` 分配 `csf_counts_by_block(0:block_count)`，第 0 项为 0；不能调用当前会重新分配并覆盖块状态的
  `install_descriptor_header`。全局 CSF 数在加载期间始终指全基组，块前缀和用局部宽整数计算。
- RMCDHF 保留全基组占据数及每块有符号 J/宇称信息的低内存布局，避免从约一份 byte 数组膨胀到五份。
- 对照 core、轨道数、电子数、总 CSF 数、逐块长度、J/宇称、rank 配置；保留 MCP 的 DIAG/ICCUT/LFORDR。
- 早期计划记录的 `getaldwt::CASE(5) → itjpo → JCUPA` 依赖仍需单独复现和核定。
  交接时两个 MPI 源码仍调用 `itjpo`，但其 loader 不长期分配耦合数组。
  此项由本迁移任务负责，必要时以独立正确性变更和权重测试处理；不为绕过问题恢复全量耦合数组。

## 4. 输入与 MCP/混合系数身份

header 哈希只绑定 header，不绑定描述符 payload 或 CSF 顺序。
现有 MCP 也没有完整的逐 CSF 身份摘要；计数、对称性相同不能证明对应同一基组。

近期接入保持 MCP 格式：强制检查已有可比元数据；回归使用同一固定描述符生成 MCP，
记录描述符、文本对照、MCP 文件集、程序和交互输入的 SHA-256，并按 rank 数生成对应分片。
缺少可靠来源的旧 MCP 应重新生成；不能给历史文件补一个文件名就宣称完成配对。

完整 payload/顺序摘要与 MCP、mixing 产物的绑定列为单独协议工作，核定版本、旧文件兼容和消费者后实施。
不把尚未写入产物的身份字段作为当前 loader 能验证的能力。

## 5. 分阶段推进与生产入口退役

| 阶段 | 工作及交付 |
| --- | --- |
| D0 盘点和基线 | 完善 §1.3 的所有消费者/产物/重读路径表，记录原输入与数值数组需求；冻结合法文本、V2 Parquet/header 和旧 GCSDESC1 对照夹具、基准程序及 SHA-256；确认 Rust 依赖和模块隔离方案 |
| D1a Rust reader 与 C ABI | 在 rCSFs 提取无输出副作用的分批 reader、严格 header/metadata 校验、元数据预扫描及 header-only；新增可选无 Python FFI，定稿 ABI/所有权/错误和静态库构建；含纯闭核、失败重开及双句柄测试 |
| D1b Fortran 公共适配 | 链接 Rust，替换二进制文件解码层，打通深层错误返回、三数组和仅占据数接口；接入并复验 rangular 串行/MPI 的直接 Parquet 输入，不再以旧 GCSDESC1 入口通过测试替代 |
| D2 RMCDHF | 优先完成两个 MPI 求解器，再覆盖两个串行版本；共享解码/MCP 对照，保持低内存和各自初始化契约；单独验证标准权重路径 |
| D3 其余应用 | RCI、biotransform、transition/phase、JJ2LSJ、密度/性质及 rwfnestimate 全部接入；重点验证双基组、轨道映射、变换后输入和程序之间的描述符产物链 |
| D4 唯一生产输入契约 | 为全部程序确定 Parquet/header 的默认命名、显式选项和双基组参数；移除文本选项/交互提示、自动探测与回退，旧格式明确非零拒绝；同步 CLI、脚本、工作流和使用说明 |
| D5 删除旧实现与最终验收 | 核对无生产调用后删除 suite 旧 CSF reader 和 GCSDESC1 reader；退出生产使用的 Python 预处理脚本移除或明确归档为离线历史工具；清理 CMake/安装产物/测试中的隐性依赖，逐程序完成无 Python、无旧文本、无 GCSDESC1 验收 |

**全部程序只读描述符是本任务完成条件，不是可延后的愿景。** 可选描述符入口或只迁移 RMCDHF 都不能
宣布全任务完成。阶段性并存仅用于实现过程和对照；不能在仍有应用依赖时提前删除共享文本 reader，
也不能在最终发布中保留 `--csf-text`、扩展名/内容自动回退或恢复文本的隐藏路径。

实施更新：D0 已冻结独立 HEAD 基线；D1a/D1b 的 reader/ABI/公共适配、D2 四个 RMCDHF、
D3 其余消费者、D4 默认/显式入口、D5 生产 reader 退役均已落地。
Check/Release/无 MPI/MPI-only 的功能矩阵已通过；不将其写成历史告警门禁或物理精度认证。
旧文件可以由独立的原生 rCSFs 工具离线转换，计算程序不自动转换；对照文本保留为测试数据或外部基准输入。
文本输出仅在确有外部互操作用途时保留并列清单，不得成为 suite 程序之间的必需输入。
不要求删除 rCSFs 的通用文本转换/恢复功能，也不把径向文件、mixing、MCP、标签文件误当作本次要退役的 CSF 文本。

旧 CSF 路径的最低限度保护归本任务：失败后不得消费记录、必须检查解析返回值、错误必须非零退出，
并保留合法闭核语义。若本批直接替换入口，将已有失败夹具转成新 reader 的验收；
若旧入口继续使用，先做最小保护，不全面重构将退役的文本解析器。
对应 P1 发现仍有效，并非因未来迁移就视为已修复。

## 6. 验证与完成条件

- Rust 格式层：直接读取真实 rCSFs 生产的 V2 Parquet/header；覆盖所支持的压缩、字典编码、多个 row group、
  footer-only 元数据、强制小批次跨块/row group、零 peel、相邻同对称性块。确认不全量物化 rows。
- 严格拒绝：缺失/重复/冲突元数据、错误 header 字节绑定、列错序/错类型/实际空值、旧 V1/printed-v/归一化数据、
  普通三行 CSF Parquet、损坏 footer/page/压缩流、截断、计数/偏移/尺寸溢出、轨道重复和非法 core/peel。
  旧 GCSDESC1 的端序/长度测试只用于过渡基线，不能替代 Parquet 格式测试。
- FFI 生命周期：C ABI 与 Fortran 驱动共同核对类型宽度/布局、缓冲区容量、批尾和 EOF、状态码、错误消息、
  初始化失败/读中途失败清理、关闭后重开、左右态两个句柄交错读取及可控分配失败；明确不可恢复错误边界。
- 公共 reader：元数据、数组逐项精确对照；空/满壳、开放 J=0、seniority、多个及相邻同对称性块、完整闭核、顺序保留。
- 物理拒绝：非法态/耦合/电子数/宇称、错误 block、无法唯一补全的 J/K、suite 不支持的额外 w 或轨道/状态范围。
  Rust 恢复成功不能替代 Fortran 数值数组语义一致性验证。
- MPI：1/2/4 ranks、空工作 rank、强制小消息、root 读取错误、worker 分配/契约错误，均明确失败且无挂起。
- RMCDHF：占据数、逐块 J/宇称、下界和 MCP 元数据一致；不增加全量耦合数组；标准权重路径另有解析期望测试。
- 应用集成：每个覆盖表条目在只提供 Parquet/header、没有旧 CSF 文本和 GCSDESC1、没有 Python 解释器/包的
  运行环境中完成；特别验证双基组、变换后重读和程序间产物链。迁移期间与固定基准程序在相同配置、rank 数
  和交互输入下比较正常输出、迭代数和停止原因；最终测试可以使用冻结的期望数组/产物，不保留生产文本 reader。
  原生 reader 的退出码 0 不能单独作为成功依据，须验证完成标志与完整输出。
- 入口拒绝：仅有 `.c`/`rcsf.inp`、错误描述符旁放合法旧文本、仅有 GCSDESC1，均不得成功或自动回退。
- 构建部署：无 Python feature 的 Rust reader/FFI、GNU 15 Release/Check、受影响的 MPI-only/serial-only 配置和
  当前 CTest 清单；核查 Python/PyO3 和 CLI 依赖未进入部署库、静态库链接需求及离线构建；
  不恢复已移除目标或放宽历史告警门槛。
- 退役审计：全源码/脚本/CMake/文档扫描加调用链核对，列出并解释剩余文本引用；只允许离线转换、测试夹具、
  明确的外部互操作输出或不属于 CSF 输入的文件。不得用“默认已切换”代替“旧入口已删除”。
- 仅在用户授权的运行范围内做本地小型软件回归；另行开始性能/科学实验前说明假设、输入、对照和判据，不提交 SLURM 作业。

旧 14/6 项描述符专项、67/77 项全套回归和历史计时是证据索引，不是本次验证结果或固定验收数量。
Python 可以继续用于开发测试和基线比较，不属于最终生产输入依赖；suite 测试使用现有共享 Python 环境，
rCSFs 的 Python 测试遵循其独立环境规范。构建使用锁定的 Rust 工具链、GNU Fortran 15，
并显式设置 `OMPI_FC=/usr/bin/gfortran-15`；shell 示例采用 fish。

完成每阶段后更新应用覆盖表、输入契约、回归证据和剩余文本调用。
现有共享库不复制，CSFG/退役求解器/rangular OpenMP 不恢复。
用户另行要求前不提交、不推送，不覆盖已安装程序。

## 7. 与其他任务的交接及启动提示

本任务涉及 rCSFs 的 reader/可选 C ABI 和 suite 的 CSF 适配模块、旧 reader 的过渡保护及 RMCDHF 输入状态依赖。
修改分别留在拥有行为的仓库；rCSFs 的 native-reader 和 suite 的共享输入模块已新增。
重命名任务暂缓这些文件的大规模改名；静态诊断任务将原 CSF P1 标为“由迁移任务处理”，不重复实施。
RIS4 角向缓存和 JJ2LSJ 唯一标签文件不是 CSF 输入，仍归静态诊断任务。

开始时先阅读两级 `AGENTS.md` 和当前 `git status`，保留所有无关改动。
多对话同时改同一仓库时按文件排程或使用独立检出，合并后重新核对来源快照和测试清单。

> 请按 `grasp-atomic-suite/docs/plans/csf_descriptor_migration_2026-11-16.md` 启动 CSF 描述符读取迁移。
> 采用 Rust 原生 Parquet/header reader，经版本化 C ABI 进入 Fortran；先按本文复用盘点在 rCSFs 提取 reader，
> 补严格绑定、纯闭核、分批生命周期、无 Python 构建及 FFI，再接入 suite 公共层和 rangular 串行/MPI。
> 随后优先完成普通/memory RMCDHF MPI，再覆盖两个串行版本及所有其余消费者。
> 最终只保留描述符生产入口，删除旧 CSF 文本和 GCSDESC1 输入依赖，不调用 Python/CLI 转换或恢复临时文本。
> 保留块顺序、低内存布局、MCP 契约、双基组映射和原数值流程；旧 MPI 计划只用于仍适用的细节，不恢复退役目标。
> 同步处理过渡期仍在用的旧 CSF 错误路径并完成必要回归，逐阶段更新覆盖和部署证据。
> 保留其他任务改动，暂时不要提交或推送。
