# rangular_mpi MCP 存储与无损缩减分析

日期：2026-11-16。状态：**源码核查和已有文件压缩实验完成；未修改生产格式或算法，未运行新的 GRASP/MPI/SCF 计算。**

本文的多通道 MPI MCP 布局是历史基线。Suite `0f35673` 实施了带目录的 MPI 格式。
原压缩测试保留其样本和版本范围。[现行 MCP 说明](../implemented/design-0020-rmcdhf-mcp-io.md)

源码基准：本仓库 `caaf3a51ce12ec8e85e30622fc854c9e2a070683`，工作区上游参考
`grasp` 为 `63f2e1b29975337a136bcd95eb985ffef8507e89`。范围是普通
`rangular_mpi` 及其普通 RMCDHF 消费端，不包含 CSFG 格式。
上游目录名为 `rangular90_mpi`，本仓库为 `rangular_mpi`；直接比较发现
`setmcp.f90`、`sort.f90` 相同，核心生成/排序布局相同；suite 另有计时、输入适配、
接口声明及空 rank 的 `SETSDA` 修复，不改变此处讨论的正常 MCP 布局。

## 结论

1. 文件数由 MPI rank 和角动量通道相乘：`P × (KMAX+3)`。各 J/宇称 block 追加到同一批文件，
   不再乘 block 数。减少 rank 能减少文件数，但不会按比例减少物理贡献总数。
2. 最终文件已经保存稀疏结构、只存一个三角区、按积分标签分组；不能把它视为大量零元素的稠密矩阵。
   大文件主要是显式角系数贡献的数量，而非文件头或 Fortran 记录标记。
3. 无损压缩确实有效：本轮已有串行样本 635.79 MiB，Zstd 3 级后 257.53 MiB，
   9 级后 226.00 MiB；所有解压流 SHA-256 与原件一致。此结果不是生产 MPI 大模型压缩率承诺。
4. 新格式可以去掉可重建的 `ICLMN`：默认类型下贡献 payload 从 16 降到 12 字节，
   精确保留索引、浮点位模式和顺序；本样本全套理论减少 23.08%。此方案尚未实现。
5. 若要求计算时也少占空间，优先设计“每 rank 一个有目录的容器 + 分块无损压缩 + 共享读取接口”。
   单纯归档只减少保存/传输空间，旧程序仍需解压；文件系统透明压缩可作为基础设施允许时的替代路线。

## 1. 存储的物理对象及分工

MCP 保存的是构造 Dirac–Coulomb 矩阵所需的**角系数及其稀疏位置**，不是径向轨道、
已算好的径向积分、ASF 混合系数或最终矩阵数值。概念上，矩阵元由一电子积分乘 T 系数、
二电子 Slater 积分乘 V(k) 系数，加上从占据等信息直接构造的对角贡献得到。
因此同一 CSF 对可以对应多个积分标签和多个 k，其存储量通常远大于一个矩阵元。

| 文件（rank r） | 内容 | 不包含的内容 |
| --- | --- | --- |
| `mcprrr.30` | 全局元数据、各 block 的稀疏端点 `IENDC` 和行号 `IROW` | 矩阵元浮点值 |
| `mcprrr.31` | 显式一电子 T 角系数，按轨道对标签分组 | 径向积分 I(a,b) 的数值 |
| `mcprrr.(32+k)` | 显式二电子 V(k) 角系数，按四轨道标签分组，`k=0…KMAX` | Slater 径向积分数值 |

依据：[SETMCP](../../src/apps/rangular_mpi/connect_mcp_files.f90)，
[MCPMPI](../../src/apps/rangular_mpi/compute_mcp_angular_coefficients_mpi.f90)，
[SETHAM](../../src/apps/rmcdhf_mpi/assemble_hamiltonian_block.f90)。

生成端按 `JA=rank+1, rank+1+P, …` 分配 CSF 列，内层只遍历 `JB=1…JA`。
`npos` 是当前 rank、当前 block 的稀疏位置编号；每列保留对角位置，非对角位置只在有保存贡献时出现。
不同 rank 保存不同列，通常不是 P 份重复的完整矩阵。

已有缩减包括：`ABS(coefficient)>1e-10` 筛选；一电子仅显式写满足条件的非对角 T；
排除 `k=0, a=c, b=d` 的 F0 项；四轨道标签利用积分对称性规范化。
`SETHAM` 在读取 MCP 前用占据、`FCO/GCO` 等构造无需 MCP 的对角项。
因此不能再把这些已隐式处理的项当成新压缩收益。
上述 cutoff 是当前实现事实，不是本报告推荐的精度阈值。

`KMAX=max(NKJ)`，其中 `NKJ=2j`；实际创建 `.30…(32+KMAX)`，故每 rank 有 `KMAX+3` 个文件。
例如 `KMAX=9` 时，8/64/128 ranks 分别有 96/768/1536 个最终文件。
无贡献的 k 也写合法空 block；不能直接删除。最高 k 是否为空取决于实际轨道集合，不能一律减掉 KMAX。
源码依据：[SETMCP](../../src/apps/rangular_mpi/connect_mcp_files.f90)、
[轨道量子数解析](https://github.com/YenochQin/grasp-atomic-suite/blob/caaf3a51ce12ec8e85e30622fc854c9e2a070683/src/lib/lib9290/prsrsl.f90)、
[block 主循环](../../src/apps/rangular_mpi/rangular_mpi_main.f90)。

## 2. 二进制布局

这些文件通过 `FORM='UNFORMATTED'`、默认 sequential access 写出；不是文本，也不是无边界的裸数组。
通常一个逻辑记录为 `[长度标记][payload][长度标记]`。GNU 默认标记各 4 字节；
单个超大逻辑记录可能拆为多个子记录，不能简单用“文件总长是 4 的倍数”判断有效性。
端序、INTEGER/LOGICAL 宽度和记录标记受编译器及选项影响；本轮样本为小端、4 字节整数/逻辑量、
4 字节标记、8 字节浮点。[GNU Fortran 官方格式](https://gcc.gnu.org/onlinedocs/gfortran/File-format-of-unformatted-sequential-files.html)
及 [OPENFL](../../src/lib/core/file_open.f90) 是依据。

下列每一行表示一次逻辑记录，而非每个括号成员各一个记录。

`.30` 的前置记录：

```text
NCORE, NBLOCK, KMAX
NCFBLK(1:NBLOCK)
IDBLK(1:NBLOCK)                   # 每个标识固定 8 字节
```

所有文件随后都有三条公共记录：

```text
'MCP', NBLOCK, MYID, NPROCS       # 'MCP' 为 3 字节字符
NELEC, NCF_total, NW
DIAG, ICCUT(1), LFORDR
```

注意公共头只记录第一个 block 的 `ICCUT`；它不是全部交互设置的完整 provenance。
`KEYORB` 也未写入文件；当前标签基数为 215。源码：
[SETMCP](../../src/apps/rangular_mpi/connect_mcp_files.f90)、
[SETMCPMPI](../../src/apps/rangular_mpi/open_mcp_outputs_mpi.f90)。

`.30` 随后为每个 block 追加：

```text
'MCP', block_id, NCF_block
NNONZ_local
IENDC(本 rank 的列), IROW(1:NNONZ_local)
```

`IENDC` 是累计终点，不是每列长度；列编号由 rank 和 P 隐含给出。
`IROW(index)` 给出另一 CSF 编号。默认 4 字节整数下，主要大小为
`4 × (local_columns + NNONZ_local)`。源码：[SETSDA](../../src/apps/rangular_mpi/write_mcp_sparsity.f90)。

`.31` 和 `.32+k` 为每个 block 追加：

```text
'MCP', block_id, NCF_block, NCOEFF
LABEL, NCONTR
(ICLMN, INDEX, COEFF) × NCONTR    # 同一 LABEL 的所有贡献共用一个记录
...                             # 下一积分标签
0, 0                            # block 结束；空 block 也保留
```

`ICLMN=JA`，`INDEX=npos`，`COEFF` 为 double。T 的标签是 `min(a,b)*KEYORB+max(a,b)`；
V 的标签把规范化后的四轨道号以 KEYORB 为基数打包。k 在文件编号中，不再逐项保存。
排序仅保证 LABEL 递增，**不保证同一 LABEL 中 INDEX 有序**。
源码：[MCPMPI](../../src/apps/rangular_mpi/compute_mcp_angular_coefficients_mpi.f90)、
[SORT](../../src/apps/rangular_mpi/sort_mcp_coefficients.f90)。

以 4 字节整数/标记为例，n 个贡献的 payload 是 `16n` 字节；
每个非空积分组另有 24 字节（标签/计数记录 16，贡献记录两端标记 8）。
这解释了为何最终文件不能按“每系数一对记录标记”估算。

## 3. 临时文件与峰值空间

流程是 `SETTMP → 逐 CSF 对写 tmp → SETSDA/SORT → 追加 mcp → 删除当前 tmp`。
文件名 `tmprrr.30…` 与最终 MCP 分开，但 Fortran unit 编号仍用 30、31、32+k；
不能仅凭 `WRITE(31)` 认定写入最终 `.31`。

生成阶段，每个临时系数记录为 `(JA,npos,LABEL,COEFF)`：payload 20 字节，加标记共 **28 字节/贡献**。
临时结构记录 `(JA,JB,npos)` 共 **20 字节/位置**，最终结构主要只需 4 字节/位置和列端点。
排序过程中，当前 block 的临时文件会与已生成的最终文件共存；峰值取决于 block 分布和排序进度，
不能简单等同最终文件总量，也不能把“最终文件压缩”当成解决所有生成阶段空间问题。

`SORT` 先把当前通道的全部系数、列号、位置、标签读入内存；k=1 还分配额外整数数组。
因此纯粹把临时文件全换成内存会把磁盘压力转移到 RAM。
批量写临时记录可将主数据从 28 降到接近 20 字节/贡献，约少 28.6%，还减少 I/O 调用；
但这是修改临时读写接口后的理论收益，尚未实现或测量。
大模型可考虑有界缓冲与分段外排序；只改最终输出容器不会消除此瓶颈。
依据：[SETTMP](../../src/apps/rangular_mpi/open_temporary_mcp_channels.f90)、
[MCPMPI](../../src/apps/rangular_mpi/compute_mcp_angular_coefficients_mpi.f90)、[SORT](../../src/apps/rangular_mpi/sort_mcp_coefficients.f90)。

## 4. 本轮已有文件测量

实验目的：确认当前真实字节布局的主要占比，并验证通用无损压缩是否有实际收益。
假设：整数索引及角系数位模式存在可压缩重复；收益大小事先未知。
判据：压缩后总字节数下降，且所有文件解压 SHA-256 与原件逐一一致。
未设科学误差阈值，因为本实验不改任何数值。

输入为本地已保留的一套串行输出：28 电子、41 轨道、23,459 CSFs、5 个 block、KMAX=9、12 文件。
只检查了文件及其头部，不将该样本当作任何原子模型的物理有效性证据。
没有找到本轮可用的完整 MPI 文件集；MPI 与串行的主体记录布局相同，但 rank 切分、标签分布和压缩上下文不同。

| 部分 | 字节数 | 总量占比 |
| --- | ---: | ---: |
| `.30` 稀疏结构及其头 | 43,843,538 | 6.58% |
| 38,461,197 条系数的 16 字节 payload | 615,379,152 | 92.31% |
| 系数文件的标签、计数、文件头、记录标记等 | 7,447,374 | 1.12% |
| 全部 | 666,670,064 | 100% |

另有 10,937,337 个稀疏位置、310,189 个非空积分组。`.41` 是全空通道，只有 258 字节；
删除它几乎不节省容量，却破坏旧读取契约。对全部贡献逐组检查，同一 `(block,k,LABEL,CSF对)`
没有重复位置；“合并重复贡献”对这套样本没有收益，不能假设所有大文件都含大量可合并条目。

Zstd CLI 1.5.7，每文件独立压缩、单线程 `-T1`，三个等级各测一遍：

| 方式 | 压缩后字节 | 压缩后 MiB | 减少 | 压缩秒 | 解压并 SHA-256 秒 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 原件 | 666,670,064 | 635.79 | — | — | — |
| Zstd 1 | 335,932,850 | 320.37 | 49.61% | 1.42 | 0.59 |
| Zstd 3 | 270,035,743 | 257.53 | 59.49% | 2.44 | 0.67 |
| Zstd 9 | 236,979,461 | 226.00 | 64.45% | 6.83 | 0.64 |

全部 36 次文件 roundtrip 的 SHA-256 一致。临时压缩副本校验后删除，原始文件未改。
时间为单机缓存/临时 tmpfs 条件，包含进程启动及校验，没有 fsync；不代表冷盘、共享文件系统、
实际 RMCDHF 重复扫描或端到端计算加速。未做 gzip/xz/LZ4 优劣排名。
Zstd 为无损算法；等级是压缩速度/体积的取舍，并非保留浮点位数。
[官方说明](https://meta.github.io/zstd/)及[格式规范](https://github.com/facebook/zstd/blob/dev/doc/zstd_compression_format.md)。

本轮交付附带 `mcp_storage_study.py/.json/.log`：JSON 保留每文件 hash、尺寸、计数和测量结果，
脚本记录本地输入路径与操作；这些本机材料不纳入源码仓库。不能把本表的百分比直接套用于其他模型。

## 5. 读取端决定了压缩接入方式

普通 `rmcdhf_mpi` 每次装配矩阵重绕 MCP；为每个待优化轨道构造系数时又反复扫描，
并按标签跳过不需要的贡献记录，二电子各 k 还会重绕 `.30`。
依据：[MATRIXMPI](../../src/apps/rmcdhf_mpi/solve_hamiltonian_blocks_mpi.f90)、
[SETCOF](../../src/apps/rmcdhf_mpi/assemble_orbital_potential_coefficients.f90)。因此：

- `.zst/.gz` 改名为 MCP 不能读取；旧 Fortran READ 没有解码器。
- 普通解压管道不能无修改支持 `REWIND`；每次从头解压还可能抵消 I/O 收益。
- 新读取层应支持逻辑流重置、block 定位、标签与 payload 读取及跳过；压缩块须有偏移/长度，
  才能跳过无关数据而无需全流解码。
- `.30` 结构值得缓存，系数适合有界分块解码；不应强迫所有系数常驻内存。

当前 `SETMCP` 强制要求文件头 `MYID/NPROCS` 与运行相同，结构也按同样列分工读取。
所以不能直接改头、改文件名、拼接文件来改变 rank 数。
[检查位置](../../src/apps/rmcdhf_mpi/connect_mcp_files.f90)。

`rmcdhf_mem_mpi` 启动时通过 `SETMCP_MEM` 两遍读取后缓存数据，适合先验证启动解码接口；
但现有内存版本身不减少磁盘 MCP，且 RAM 要容纳所有角系数。
它与普通版还有其他代码差异，不能把切换求解器当作已经验证的纯存储优化。
依据：[SETMCP_MEM](../../src/apps/rmcdhf_mem_mpi/cache_mcp_files.f90)、
[现有四版本比较](evidence-0007-rmcdhf-memory-variant-comparison.md)。

普通 RCI/RCI_MPI 直接调用 `ONESCALAR/RKCO_GG` 构造角系数，不消费这套 rangular MCP；
此次格式适配重点是四个普通 RMCDHF 版本，不必据此修改 RCI `.res`。
依据：[RCI SETHAM](../../src/apps/rci/assemble_hamiltonian_block.f90)、
[RCI_MPI SETHAM](../../src/apps/rci_mpi/assemble_hamiltonian_block.f90)。

## 6. 可选措施及边界

| 措施 | 文件数 | 实际数据量 | 旧程序直接使用 | 判断 |
| --- | --- | --- | --- | --- |
| tar/每 rank 打包 | 归档时减少 | 不压缩则基本不变 | 需解包 | 解决目录/传输管理 |
| Zstd 无损归档 | 单独压缩不减数量；配合打包可减 | 本样本明显下降 | 需解压 | 可用于保存；不降低恢复后的运行峰值 |
| 文件系统透明压缩 | 不变 | 物理占用可下降 | 保持原文件接口 | 取决于服务器文件系统支持，需实测 |
| 每 rank 一容器，不压缩 | `P(KMAX+3) → P` | 基本不变 | 需读写接口 | 先解决文件数及元数据操作 |
| 容器 + 分块 Zstd | 可降为 P | 有望下降 | 需读写接口 | 推荐的持续运行方案，尚未实现 |
| 删除可重建 ICLMN | 不决定文件数 | 贡献 payload 少 25% | 需读写接口 | 精确结构缩减，需测恢复成本 |
| float32/更大 cutoff/删除小项 | 不决定文件数 | 可下降 | 格式/数值契约改变 | 有损近似，不作为首选压缩 |

透明压缩例如 Btrfs 的 extent 压缩，应用看到的逻辑字节与文件名保持不变；
但不减少 inode 数，已有文件不会因启用设置自动全部改写，压缩率及吞吐须在目标存储上验证。
本轮未修改文件系统配置，也未验证集群支持。
[Btrfs 官方说明](https://btrfs.readthedocs.io/en/latest/Compression.html)。

### 精确删除 ICLMN

给定 `INDEX`，从同一 rank、同一 block 且通过结构检查的合法 `.30` 的累计端点找到所属列，
即可恢复 `JA=ICLMN`；`IROW(INDEX)` 恢复 JB。
`SETHAM` 本来把列号读入 dummy 后丢弃，但 `SETCOF` 会使用列号，因此必须提供恢复逻辑。
现有 [MCP 解析器](../../scripts/mcp_compare.py) 已用端点查找还原列号并核对该字段。

保持整数默认宽度，`(int32,int32,float64)` 可变为 `(int32,float64)`：
本样本理论节省 `4×38,461,197=153,844,788` 字节，整套从 635.79 降至约 489.07 MiB。
这是未压缩格式的算术估计，不能与 Zstd 实测收益直接相加或相乘。
恢复可选择端点二分，或额外 `index→column` 映射，需权衡 CPU 与 RAM。

### 后续可测的编码

- 把交错整数/浮点变为分字段数组，做可逆 byte shuffle，再压缩；数据顺序和 float64 位模式保持一致。
- `.30` 的行号在每列内递增，端点也递增，适合差分/变长整数；系数 INDEX 则不能假定有序。
  如对系数用差分，应保持原序并支持有符号差值，或显式保存恢复顺序。
- 对逐位相同的浮点值做字典编码，在组内有 n 项、u 个不同位模式、索引宽度 b 字节时，
  仅当 `8u+n*b+元数据 < 8n` 才启用。相同值在不同 CSF 对上仍是不同贡献，不能删除。
- 去掉全部记录标记并不会消除主体贡献；本样本系数文件所有非贡献开销合计仅 1.12%，
  不值得把“去 Fortran 标记”作为主要容量方案。

上述均未做生产实现。重新排序或合并贡献可能改变浮点累加顺序；
按容差合并浮点值、降低精度、提高 cutoff、DIAG/LFORDR 截断会改变数值/模型，应单独论证。
新格式还应保留大尺寸能力：不能为了缩小文件把可能溢出的稀疏位置和计数强制压到 16 位；
大容器偏移应明确使用 64 位。现有默认整数计数的溢出风险也不会因压缩自动消失。

## 7. 建议实施顺序及验收目标

1. **首先取得生产 MPI 样本的完整存储分布**：每 rank、k、block 的贡献数、NNONZ、字节、
   临时峰值及普通 RMCDHF 的读取时间；用同样逐文件 roundtrip 方法验证压缩率。
2. **建立共用 MCP I/O 接口并保留旧格式**，先保持所有字段和贡献顺序，
   用每 rank 一个容器消除通道文件的乘数。独立 rank 文件可避免首版就引入全局共享写锁或 MPI-IO 集中竞争。
3. **加分块无损压缩**：目录记录 `(rank,channel,block,label/segment)`、偏移、压缩/原始长度、计数、校验；
   大标签也拆成有界块，并支持无收益时原样存储。显式 magic/version、端序、整数宽度、float64、KEYORB、
   CSF/轨道顺序身份及分区元数据。不要每一条系数单独建压缩帧。
4. **再加精确结构缩减和临时 I/O 批量化**，分别测量收益；每阶段能还原旧格式、回退并定位差异。
   跨进程数重分区是另一个功能，需要重建每 rank `.30` 和局部 INDEX，不随容器合并自动完成。

验收先比较还原后的完整字节/字段、稀疏结构、原始贡献顺序和 float64 位模式；
覆盖多 block、空 rank、空 k、大记录、截断/校验失败、重绕与跳过。
已有 [mcp_compare](../implemented/design-0024-mcp-compare.md) 可做结构和零容差贡献比较，但它会对齐贡献，
不能单独证明文件字节或原始累加顺序不变；仍需独立 roundtrip/hash 检查。
然后用相同输入、初始轨道和 MPI 数对照 RMCDHF 矩阵、SCF 轨道/能量及实际 I/O、RSS、峰值磁盘空间。

本轮支持“存在无损容量收益及可重建字段”的判断；未验证新容器、服务器吞吐、
MPI 分块压缩或任何新的科学计算结果。下一步是生产 MPI 样本测量与存储接口原型，
不需要先改变角系数精度或物理模型。
