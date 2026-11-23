# RMCDHF MCP 阶段 3 实施过程

记录日期为 2026-11-23。A 的源码基准是 Suite `c1377aa0`，C 的源码基准是 `0f356733`。本次交付包含独立 A/C
测试包、矩阵读盘分段及逐轮页驻留诊断。已经完成服务器和本机的数值与墙钟对照，目标盘上的矩阵分段、实际读盘
减量、B 路径和资源预算仍待验收。[对应计划](plan-0005-rmcdhf-mcp-production-io.md) · [实测记录](../reviews/evidence-0025-rmcdhf-mcp-phase3-benchmark.md)

## 实施顺序

1. **取得服务器 A/C 基线。** 作业 14 在目标 XFS 上用单节点 48 ranks 跑 Zn I AS2。冷、热态各三组，均到三轮
   EOL 上限。C 的全 SCF 最慢 rank 墙钟缩短约 10.7%，冷态全 SCF `read_bytes` 增加 94.32 MiB。作业 15 在两套
   构建的相同调用位置加入势系数窗口和 MPI 集合调用计时。势系数墙钟缩短约 94.4%，该窗口两版均测得零读盘，
   全 SCF 的额外读量仍待按矩阵及其他窗口定位。
2. **确定配套构建方式。** 搜查了本地 Git 历史、交接归档和服务器测试目录，只找到可核验的 A/C 构建。按计划使用
   A/C 备用方案，两版分别生成 MCP。矩阵分段包改为直接导入已有的 Suite `libcsf_descriptor.a`。构建脚本在缺少库时
   报出路径，摘要记录实际链接的库哈希。两份历史源码单独导出，采样代码只注入基准快照。
3. **补齐矩阵窗口和本机重复。** A 的采样包含旧文件头重定位、稀疏结构读取和未选根块的顺序消费，C 的采样涵盖
   带目录结构和矩阵系数读取。本机用 8 ranks 在同一 NVMe 的 Btrfs、XFS 分区各做三组冷、热配对。C 在两种配置下
   的全 SCF 均缩短约 44%–46%。冷态读量随配对变化。生成期另采样文件逻辑大小和已分配空间，同一版本在两种文件系统上
   的最终 MCP 逻辑大小一致。XFS 一次热态预检查失败，确认顺序预热后早期分片已有页被回收。保留已完成的运行和失败摘要，
   对剩余热态改用 `mmap` 逐页预热，验证全部驻留后续跑。
4. **从首次采样扩展到逐轮定位。** 首版诊断记录初始矩阵装配和首次势系数构造前后的 `mincore` 页驻留和 `read_bytes`。
   本机 XFS 首次诊断中两版完成初始矩阵后均全量驻留，势系数读盘均为零。逐轮版按 EOL 矩阵求解编号，每轮记录
   首次构造前后、首次正读量的调用序号和该轮累计读量。用实际冷读 16 MiB 的受控测试验证了第 2 次调用的事件定位。
5. **验证逐轮诊断并归档。** 本机 XFS 的 8-rank 冷态 A/C 复跑中，C 在首轮矩阵末块装配后驻留 17.74%，
   至首次势系数构造前降至 4.76%。八个 ranks 均在第一轮第 1 次调用发生读盘，共 3,857.53 MiB，后两轮均为零。
   A 的三轮势系数窗口均为零。两版迭代能量、停止状态和输出哈希与同 rank 数的非诊断对照一致。该采样定位了读盘
   窗口，页回收来源仍待结合设备与内存计数核查。

## 构建与复跑入口

| 用途 | 交接包 | 运行与验收 |
| --- | --- | --- |
| 矩阵分段性能对照 | `stage3-matrix-external-kit.tar.gz` | `build.fish` 导入外部 Rust 静态库。按 `SLURM.md` 在目标盘复跑冷、热配对，核对 `summary.json` 中的矩阵、势系数和其余 SCF 窗口，再比较逐组墙钟。 |
| 逐轮页驻留诊断 | `stage3-cache-epoch-external-kit.tar.gz` | 本机用 `run-local.fish` 指定工作盘和归档目录，服务器按 `SLURM.md` 分别提交冷态和热态。对照 `mcp_cache_samples.epochs` 及原始 `mcp-cache.csv`，核查每轮读量与 `scf-perf.csv` 一致。 |

两份包都附 SHA-256 文件，使用原样的六份 Zn I AS2 输入。两套 Release 构建经 GNU Fortran 16.2.1 本机复核。
逐轮诊断包的双 rank 小样本覆盖冷、热和关闭探针三种情况，`.mix`、`.wfn` 与原样本一致。正式性能计时使用矩阵分段包，
页驻留诊断单独运行。探针的 `timestamp_ns` 可与外部采样的 UTC 对齐，进程累计读量和设备累计读量均按采样窗口取差。

## 证据位置

下列路径以父工作区为基准。本机工作盘保留 MCP 和中间文件，`temp` 归档保留摘要、计数 CSV、日志及校验后的数值输出。

| 证据 | 工作区归档位置 |
| --- | --- |
| 服务器作业 14、15 的最终摘要及交接包 | `temp/suite-dev-test/rmcdhf-mcp-phase3-handoff/` 下的 `summary-14.json` 和 `summary-potential-15-server.json` |
| 本机 Btrfs 三组冷、热配对 | `temp/suite-dev-test/rmcdhf-mcp-local/mcp-repeat-20261123-201558/paired-summary.json`，内含首组和补测摘要哈希 |
| 本机同盘 XFS/Btrfs 对照 | `temp/suite-dev-test/rmcdhf-mcp-local/mcp-xfs-20261123-212053/filesystem-comparison.json`，同目录保留缓存失败和续跑记录 |
| 本机首次页驻留诊断 | `temp/suite-dev-test/rmcdhf-mcp-local/mcp-cache-20261123-221704-377174/summary.json` |
| 本机逐轮页驻留诊断 | `temp/suite-dev-test/rmcdhf-mcp-local/mcp-cache-epoch-20261123-223449-392471/summary.json` |

## 下一次交接

先在目标 XFS 以 48 ranks 执行矩阵分段基准，确认作业 15 的额外读量出现在哪个窗口。再用逐轮诊断包和已有 MPI
外部采样器同时记录进程 I/O、设备读量、页缓存回收和内存压力，优先核对最后一个初始矩阵块到首次势系数构造返回的
时间窗。目标盘上若仍缺少 A 的势系数读盘，按计划选择实际发生读盘的生产算例。B 构建和目标作业的单 rank 内存预算
仍需补齐。[计划中的验收判据](plan-0005-rmcdhf-mcp-production-io.md#待实施与验收)
