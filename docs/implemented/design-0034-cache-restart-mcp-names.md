# 角缓存、RCI 重启和 MCP 文件命名

依据 [Suite 文件依赖核查](../reviews/evidence-0033-output-file-formats.md#后续迁移顺序) 的第 3 项实施。
下表中 `S` 是状态名，`A`、`B` 是按输入顺序排列的初末态。MPI 文件仍由各 rank 在原工作目录中访问。

| 用途 | 原文件名 | 现行文件名 | 复用和恢复条件 |
| --- | --- | --- | --- |
| 双正交一体角系数 | `S.TB` | `S.biorth.angular.cache` | 校验该侧的描述符摘要、基组尺寸、κ 分组、J 筛选和 rank 布局。 |
| 辐射跃迁角系数 | `A.B.±nT` | `L.A.B.rtransition.eN.angular.cache` 或 `L.A.B.rtransition.mN.angular.cache` | 保留原路径中的初态名字节数 `L`，例如 `4.left.right.rtransition.e1.angular.cache`。校验有序的基组摘要、电或磁多极算符、J 筛选、现有块头及 rank 布局。 |
| 同位素位移角系数 | `S.IOB`、`S.ITB` | `S.isotope.one-body.angular.cache`、`S.isotope.two-body.angular.cache` | 校验基组摘要、一体或二体算符、轨道标签基数和筛选阈值。串行和 MPI 共用完整缓存。 |
| jj–LSJ 变换 | `S.lsj.T` | `S.lsj.transform.cache` | 顺序回读，校验 V3 描述符 SHA-256、d 壳层约定、块元数据、每块实际保留的 jj-CSF 行和完成记录。重新生成旧缓存。 |
| 串行 RCI 重启 | `rci.res` | `S.rci.restart` | 核对 `R92RS2` 的基组计数和 V3 描述符 SHA-256 后，恢复核模型、Hamiltonian 选项、网格和轨道、稀疏矩阵及 Breit 位移。 |
| MPI RCI 重启 | `rciRRR.res` | `S.rci.rankRRR.restart` | 每份 `R92RS2` 分片均核对基组身份，并要求相同的 rank 数和分片布局。rank 编号至少三位。 |
| 串行 MCP | `mcp.30`、`mcp.31`、`mcp.(32+k)` | `S.mcp.30`、`S.mcp.31`、`S.mcp.(32+k)` | `rangular` 写出稀疏结构及 T/V 通道，`rmcdhf` 和 `rmcdhf_mem` 读取同一状态的整组顺序无格式文件。 |
| MPI MCP | `mcp.mcpidx` | `S.mpi.mcpidx` | `rangular_mpi` 发布的 `GMCPIDX1` 分片继续用基组摘要、写出批次、rank 布局和 `MCPDONE1` 校验。两种 MPI 求解器使用同一个命名契约。 |

旧角缓存缺席时，双正交、跃迁和 RIS4 走原有重算路径。基于相同记录格式、匹配输入和运行条件的缓存，
也可按表复制到新路径，让现有校验决定能否复用。旧跃迁名中，符号是算符宇称，
电多极为 `(-1)^N`，磁多极为 `(-1)^(N+1)`。

旧 `R92RES` 文件缺少基组摘要，Suite 会拒绝续算。请用原有输入冷启动 RCI 生成 `R92RS2` 文件。
新格式将 64 字符摘要追加到第二条顺序无格式记录的四个计数字段之后。续算要求原 V3 描述符、
相同的顺序记录 ABI 和完整的 MPI 分片集合。[RCI 续算记录](../../src/lib/core/rci_restart_records.f90)

旧 MCP 文件可在确认原计算状态及运行条件后复制到新名。串行 MCP 迁移要求原始 CSF 顺序及
相同的顺序记录 ABI，MPI MCP 要复制整套 rank 文件。旧 MPI 顺序通道应由当前 `rangular_mpi` 重新生成。
生产程序按状态名打开文件。MCP 比较工具仍可读取旧归档，并允许用 `--state` 和
`--candidate-state` 从各自目录选取状态。 [MCP 比较工具](design-0024-mcp-compare.md)

JJ2LSJ 的 `LSJTRN02` 文件先写基组身份和块头，再按计算顺序写块编号、保留行掩码及原有
`(jj 行, LS 行, 系数)` 记录，末尾写 `LSJDONE1`。回读要求相同的块序列和保留行。
变换系数可用于同一基组中满足这些条件的混合向量。带目录的随机读取仍由 MPI MCP 提供。

`rmcdhf_mpi --mcp-trace` 和 `--mcp-verify` 分别写出 `S.mcp.io.csv`、`S.mcp.bits`。
[MPI MCP 格式和诊断](design-0020-rmcdhf-mcp-io.md)
