# RMCDHF MCP 阶段 3 实测与续跑

2026-11-23。归档基线 A 为 `c1377aa0`，dev 路径 C 对应 `0f356733`。Zn I AS2 在目标
XFS/NVMe 上的第一轮生产对照已经完成。dev 的最慢 rank 在冷、热缓存下均加快约 10.7%，
冷态全 SCF 实际读盘量增加 2.37%。阶段 3 仍需势系数分段测量和单 rank 内存预算核对。
参见[阶段 3 验收定义](../plans/plan-0005-rmcdhf-mcp-production-io.md)。

## 两套基准的代码边界

- 路径 A 保留按 rank、按通道写出的旧 MCP 和逐轨道扫描。路径 C 每 rank 写一个带块和组目录的
  `mcp.mcpidx`，有选根的块按目录定位。普通 MPI 版在每次 EOL 密度权重更新后重建收缩表，
  轨道和乘子复用组收缩。MPI 内存版读取同一格式。dev 的计算入口使用目录格式单路径。
  [实现和阶段 0–2 验证](../implemented/design-0020-rmcdhf-mcp-io.md)
- 代表性多块算例对照采用同一份 V3 描述符和输入，由各自的 `rangular_mpi` 重新生成 MCP，
  各自的 `rmcdhf_mpi` 求解。阶段 3 测量代码在两个独立的源码快照中编译。
  快照生成器以 Git 提交为输入，记录提交、探针和注入脚本的 SHA-256。
- 作业 14 的探针覆盖 `run_scf_iterations_mpi` 全段。后续测试包保留全段数据，
  增加势系数构造与已有 MPI 集合调用的分段计时。两份基准都采用相同采样边界和探针源文件。

## 作业 14 的生产对照

Zn I AS2 含 4 个对称块，每块分别有 6,365、17,701、25,897、29,534 个 CSF，共 79,497 个。
同一 Linux 节点使用 48 个 ranks，Slurm 按核绑定。归档版和 dev 使用相同的描述符、初始波函数、
核数据、交互输入及 MPI/BLAS 动态库。SCF 上限为三轮。冷、热缓存各做三次配对对照，
相邻重复交换 old 和 dev 的运行顺序。每次启动 RMCDHF 前用 `mincore` 核对 MCP 文件页，
冷态均不驻留，热态均驻留。

### 三次重复的中位数

| 缓存与指标 | 归档版 A | dev C | dev 变化 |
| --- | ---: | ---: | ---: |
| 冷 · SCF 实际读盘，全部 ranks（MiB） | 3,973.14 | 4,067.46 | +2.37% |
| 冷 · 峰值 RSS，最大 rank（MiB） | 569.44 | 569.52 | +0.01% |
| 冷 · SCF 墙钟，最慢 rank（秒） | 263.17 | 234.88 | −10.75% |
| 冷 · 包含 MPI 启动的整次 RMCDHF（秒） | 265.74 | 237.46 | −10.64% |
| 热 · SCF 实际读盘，全部 ranks（MiB） | 0 | 0 | 0 MiB |
| 热 · 峰值 RSS，最大 rank（MiB） | 569.57 | 568.49 | −0.19% |
| 热 · SCF 墙钟，最慢 rank（秒） | 263.30 | 235.20 | −10.67% |
| 热 · 包含 MPI 启动的整次 RMCDHF（秒） | 265.81 | 237.78 | −10.55% |

SCF 实际读盘量是各 rank 的 Linux `/proc/self/io` 中 `read_bytes` 在 SCF 入口和返回点的差值之和。
预热和程序初始化发生在采样窗之前。峰值 RSS 由 SCF 返回后的 `VmHWM` 读取，
涵盖各 rank 从启动到采样点的内存峰值。最慢 rank 墙钟是每次运行中 rank 的最大本地 SCF 耗时。

三次冷态对照中，A 每次均读 4,166,135,808 字节，C 每次均读 4,265,037,824 字节，
每次增加 98,902,016 字节（94.32 MiB）。冷态最慢 rank 每轮节省 26.94–28.29 秒，
热态每轮节省 26.58–32.27 秒。C 在六次 SCF 中的最大单 rank RSS 为 569.64 MiB。

两版的 12 次运行具有相同的迭代号、输出能量、停止状态以及 `.mix`、`.wfn` 的 SHA-256。
各次均到达设置的三轮迭代上限。初始和三轮加权平均能量分别为
−1794.5696497、−1794.6448855、−1794.7740417、−1794.7829098 Hartree。

MCP 每版生成一次。A 用 22.62 秒生成 576 个分片、4,272,972,236 字节，
C 用 22.20 秒生成 48 个分片、4,284,016,692 字节。C 的最终大小增加 0.26%。
生成只测了一次，临时文件的磁盘峰值和重复生成成本待测。

### 验收判断

- **SCF 数值一致，墙钟可重复下降。** 冷、热态的全部六组配对都加速。
  热缓存时两版的 SCF 实际读盘均为零，仍观察到加速。逻辑扫描、收缩计算与 MPI 等待
  对加速的各自贡献需要分段计量。
- **势系数实际读盘条件待证。** 冷态全 SCF 读量增加。全 SCF 读量包含矩阵装配和势系数阶段。
  在势系数构造前，矩阵装配可能已经把相应 MCP 页读入缓存。要判断计划中的势系数物理
  读取是否下降，需要独立采样这一段。生产 `--mcp-trace` 的逻辑读取字节用于解释扫描和跳过，
  与 Linux `read_bytes` 的存储读取归因分开比较。
- **内存预算、B 对照和生成峰值待核。** 用 C 的 569.64 MiB 最大 rank 峰值与目标预算核对。
  作业 14 测了 A/C。若要拆分收缩复用与磁盘格式的增量，需在同一文件系统补测归档 B。
  串行与 MPI 内存版的性能验收各自使用匹配基线。

## 势系数补测的实现和验证

| 指标 | 探针窗口和汇总 |
| --- | --- |
| 势系数实际读盘与墙钟 | 包住每次 `assemble_orbital_potential_coefficients`。覆盖普通、追加和重试轨道更新及 Lagrange 乘子的两侧。每 rank 累加 `read_bytes` 差值和 `MPI_Wtime`。 |
| 最慢势系数 rank | 取上述累计构造墙钟的最大值，同时记录该 rank 的读量及 MPI 调用耗时。 |
| MPI 等待与通信 | 分别计径向势 `MPI_Allreduce`、乘子 `MPI_Allreduce` 和直接系数汇集的 `MPI_Gather`／`MPI_Bcast`。累计阻塞调用时间包含通信、归约和等候。集合调用沿用求解器原次序。 |
| 全 SCF 对照 | 保留全 SCF 读盘、最慢 rank 墙钟、最大 rank RSS、整次调用墙钟、缓存驻留和数值一致性检查。 |

复跑包 `stage3-potential-kit.tar.gz` 含两份带探针快照、原算例、构建及提交脚本。
每个 SCF 目录保留每 rank 的 `scf-perf.csv`，`summary.json` 保留逐轮运行与按冷、热态计算的中位数。
重点字段为 `potential_read_bytes_all_ranks`、`slowest_rank_potential_seconds`、
`potential_read_bytes_at_slowest_potential_rank`、`potential_mpi_seconds_at_slowest_potential_rank`、
`max_rank_potential_mpi_seconds` 和 `ranks[].mpi_seconds_by_stage`。
`status: passed` 表示采集完成且数值对照一致。阶段 3 性能判定仍以每轮实测差异为准。

两份 GNU Fortran 16 Release 构建和本地双 rank 小样本均通过。归档和 dev 的三类 MPI
计数一致，`.mix` 和 `.wfn` 与上一版小样本逐位一致。作业 14 的六份算例输入校验值
与补测包一致。macOS 本地验证覆盖数值与计时采集，Linux 读计数的解析用合成数据校验。
目标文件系统上的势系数补测待执行。

### Rust 静态库复用

`build.fish [old-build-lib-dir [dev-build-lib-dir]]` 为每个基准单独接收含
`libcsf_descriptor.a` 的目录。代码快照的 CMake 直接导入该静态库并沿用平台链接依赖。
提供两份兼容库后，构建只编译和链接两套 Fortran 程序。参数留空时，对应基准使用
内置 Rust 源码构建。输出摘要记录实际链接的静态库 SHA-256。

`c1377aa0` 自带的库可供 A 使用。C 需要新增的 `rcsfs_descriptor_basis_id`
来取得与 V3 描述符绑定的基组身份。目标平台的 `0f356733` 库可供 A、C 共用。
本地已经用同一份 dev 静态库完成两套构建和双 rank 数值对照。复现作业 14 时，
优先使用分别对应两版 Suite 的库。平台和体系结构应与测试节点匹配。

## 换机与继续测量

工作区 `temp/rmcdhf-mcp-phase3` 下保留作业 14 的 `summary-14.json` 和 `slurm-14.log`。
同目录下 `potential-probe` 的补测包含校验文件及服务器部署命令。换机时携带
`rmcdhf-mcp-phase3-handoff.tar.gz`，其中包含补测包、校验文件、作业 14 证据和文档补丁。
按解包后的 `README.md` 恢复文档和续跑。源代码快照由测试包提供，服务器上的 MCP 分片
与构建目录可以继续留在计算盘。

1. 在目标服务器的独立工作目录解包，执行包内 `SLURM.md` 的 fish 命令。先加载 GNU Fortran 16、
   MPI、BLAS/LAPACK 和 Python 3.11 以上版本，运行 `build.fish`，随后用单节点 48 ranks 提交。
   每次运行从 Slurm 作业 ID 取得唯一输出目录。
2. 作业在 `rangular_mpi` 生成和 SCF 阶段将详细输出写入各阶段目录。看似停滞时，
   先检查 `squeue` 的排队状态和对应阶段的 `stdout.log`、`stderr.log`。
   编译脚本已使用 `compiler_version` 以避开 fish 的只读 `version` 变量。
3. 取回新作业的 `summary.json`，检查缓存驻留和数值一致性。冷、热态分别配对分析
   势系数读量及最慢 rank 墙钟，再用 MPI 分项与全 SCF 指标解释差异。
   若冷态在势系数段两版都接近零，记录首轮矩阵装配对 MCP 页缓存的影响，再判断是否需要
   加采矩阵分项或选择另一组有重复物理读取的代表算例。
