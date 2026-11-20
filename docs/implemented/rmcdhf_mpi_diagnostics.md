# RMCDHF MPI 外部停顿采样与原版流程对照

日期：2026-11-15。源码核查：suite `3f82922d2caa4de061c7f62c1b604770d0092847`，
工作区 grasp `488e875d1a7ec631baf4ad33682798bc795b7435`。
采样器为此次工作树新增代码；随后按用户要求移除普通 `rmcdhf_mpi` 的核极化入口。

状态：Linux 外部采样器已实现，模拟 `/proc` 回归通过；未取得服务器作业，
未在真实 Linux 作业上验证采样开销、rank 覆盖或停顿根因。Fortran 改动仅移除普通 MPI
的 NP 入口及配套参数，没有修改 SCF 主循环或启动原子计算。本工具落实[诊断计划](../plans/rmcdhf_mpi_single_node_diagnostics_2026-11-09.md)
的 D0 首版；D1 内部阶段计时和 D2/D3 对照、优化仍待实施。

## 1. suite 与原版的计算流程是否一致

**普通 `rmcdhf_mpi` 保留原版主要计算流程，但不是逐字节相同的求解器。**
对应上游目录为 `src/appl/rmcdhf90_mpi/`，不要与 memory、CSFG 或已退役的 orbopt 混比。

移除 NP 入口后的逐文件比较：116 个映射文件中 99 个完全相同、17 个有差异；116 个上游文件的
hash 均与导入快照一致。运行 `test/npol/check_upstream.py`，四个普通程序合计
492 个文件在逆向撤销登记的兼容修改与 NP 入口后通过快照检查。
该检查保护来源，不是端到端数值或性能证明，也不覆盖全部共享库。

EOL 主路径仍为：

```text
MPI 启动、MCP/CSF/核与网格/初猜输入
  → 初始 MATRIXmpi → NEWCOmpi
  → 每轮 SETLAGmpi
  → 按 IORDER 顺序 IMPROVmpi（跳过冻结轨道）
  → 按自洽量进行追加改善
  → 正交化及 ORBOUT
  → MATRIXmpi → NEWCOmpi
  → 原有收敛判定或迭代上限
```

[SCFmpi](../../src/apps/rmcdhf_mpi/scfmpi.f90)、
[IMPROVmpi](../../src/apps/rmcdhf_mpi/improvmpi.f90)、
[SETLAGmpi](../../src/apps/rmcdhf_mpi/setlagmpi.f90)、SETHAM、NEWCOmpi 源文件与原版一致。
普通 MPI 的停止条件包括自洽量及相对加权能量变化，不能套用串行说明或退役变体的判据。
SETCOF 仍扫描 MCP 文件；Davidson 的 SPICMVmpi 仍按 rank 分配列并归约完整结果向量。
共享库相关路径的显式接口修改没有把它改成 OpenMP 或另一种并行算法。

必须区分以下差异：

| 差异 | 对对照的含义 |
| --- | --- |
| GETSCDmpi 核极化初始化已移除 | 2026-11-15 恢复 NUCPOT 后直接 SETRWFMPI；普通 MPI 不再读取 npol.inp，不执行 NP 检查、广播或数组复制。其他求解器的可选 NP 功能保留。 |
| MATRIXmpi 的 DDOT 改为正确双精度接口 | 原版错误声明为 INTEGER；影响触发 CI 阻尼正交化的分支，不能把修复视为纯排版。 |
| IENDC 改为 `0:NCF` 分配 | 配合实际零下标访问修复越界；原版未定义行为不应要求复现。 |
| GETALDWT 用 NCF 取代未初始化的 NCMIN | 修复 EAL 用户自定权重输入分支。 |
| SOLVE 参数 INTENT、MPI datatype handle、BLAS/LAPACK 与内部接口 | 修复参数契约、类型与有界实参；可能影响原版错误路径或编译优化。 |
| 共享 Davidson 的零列 DGEMV 边界等修复 | 消除越界实参，不改变预期收敛策略；不能据应用目录快照宣称整个库完全相同。 |

当前两个已安装二进制的编译器/MPI 链接记录也不同；源码对照不能证明它们使用同一构建配置。
若后续比较耗时，先统一 GNU Fortran 15、MPI、实际 BLAS 后端、优化选项、线程数和绑定；
记录二进制 hash、科学输入、各 block 大小、根选择、冻结轨道、网格和停止状态。
更改 rank 数必须重新生成或验证重分片后的 MCP；SETMCP 会检查 MYID/NPROCS。

## 2. 采样现有单节点 Linux 作业

入口：[sample_rmcdhf_mpi.py](../../scripts/sample_rmcdhf_mpi.py)。仅依赖 Python 标准库。
工作区使用 Tools 共享环境；独立服务器可用已有 Python 3.10+，不需要重新构建 Fortran。

在独立终端运行，从 suite 根目录执行，替换为实际本机 `mpirun` PID 和新输出目录：

```fish
../graspkit-tools/.venv/bin/python scripts/sample_rmcdhf_mpi.py \
    --pid 12345 \
    --output /path/to/local-diagnostics/run-01 \
    --interval 2 --details-interval 30 --duration 86400 --max-mib 256
```

工具只读取指定 PID 及其同用户后代，不启动、停止或修改计算，不提交 SLURM 作业。
保存 PID 与启动 tick，排除旧 PID 被重用后的无关进程；已经见过的后代在重新挂靠父进程后继续跟踪。
无法发现采样开始前已经脱离指定树的 rank，也不能保证捕获两个采样点之间短暂存在的进程。

如果启动器不在计算节点，或本机 ranks 不在同一个后代树，重复指定各本机根 PID：

```fish
python3 scripts/sample_rmcdhf_mpi.py \
    --pid 23451 --pid 23452 --output /path/to/local-diagnostics/run-02
```

不要把登录 shell、系统 init 或整台主机所有进程作为目标。多节点需要在各计算节点分别采样。
先确认 `samples.jsonl` 内的 rank 候选、exe、进程数覆盖预期；rank 环境变量缺失或冲突时保存为未知。
启动器和辅助进程也会记录，最大进程数不等于 MPI rank 数。

建议尽早附着，覆盖读取、矩阵、轨道改善和输出全过程。输出放在独立本地诊断目录，
避免向正在调查的存储设备大量写入；本工具默认不读取应用日志和科学数据内容。
输出含机器和文件路径，应保留为本地诊断数据，不提交 Git。

## 3. 记录内容与限制

| 文件/字段 | 内容 |
| --- | --- |
| `metadata.json` | 主机、内核、boot id、PID 身份、采样选项、时钟 tick 频率与页大小。 |
| `samples.jsonl` | 采集开始的 UTC、单调时钟相对时间、每次采集耗时及以下原始指标（累计计数与瞬时值）。 |
| `processes` / `tasks` | 各进程与线程状态、CPU ticks、缺页、RSS 页数、I/O 字节、上下文切换、CPU/NUMA 允许列表、上次运行 CPU。进程 CPU 总量与线程量不能重复相加。 |
| `wchan` / `stack` | 每轮线程等待位置；D 和 R 线程的内核栈按 details interval 限频。权限不足、文件消失明确保存 error，不提权。`0` 或空栈不证明没有等待。 |
| `files` | 每个进程每 30 秒记录 cwd/exe、打开文件链接与 mountinfo；不跟随链接读取文件内容。默认最多 128 个 fd。 |
| `environment` | 只保存 rank 与常见数学库线程设置白名单；不保存完整环境。设置值不证明实际线程活跃数或实际 BLAS 后端。 |
| `system` | `/proc/stat`、meminfo、vmstat、diskstats、loadavg、CPU/memory/io PSI，供时间差分析。全机计数可能含其他作业贡献。 |
| `cgroups` | cgroup v2 本组及祖先的 CPU、内存、swap、I/O 限额、事件与压力；v1 或不可见字段明确缺失。要求采样器和作业处于可对应的 mount/cgroup namespace。 |
| `summary.json` | 停止原因、采样数、出现 D 的采样数、状态观察计数、最大观测进程数和采集耗时。 |

单次快照是顺序读取，并非所有 rank 同时原子采样。`phase` 始终标记未知；
应用 stdout 有缓冲，最后一行日志也不能当作所有 rank 的实时调用栈。
缺页、I/O、CPU 与设备计数是原始累计值，分析时按时间差计算，并检查进程身份和计数重置。
`blkio_delay_ticks` 还依赖内核 delay accounting 配置，零值不能排除 I/O 等待。
PSI 的 total 单位为微秒；磁盘字段和读写等待时间按内核定义解释，不直接用 `%util` 判断瓶颈。

`--max-mib` 严格限制 `samples.jsonl`，不会写半条 JSON；另外有小型 metadata/summary 文件。
每进程最多采 1024 个线程、128 个 fd，可通过 `--max-threads` / `--max-fds` 调整；
达到限制时记录 truncated。每个 `/proc` 文本也有限长并记录截断。
达到时长/大小上限、全部跟踪进程退出或不可读、SIGINT/SIGTERM 时只结束采样器。
若使用 Ctrl-C，应在单独的采样终端中操作，避免终端同时给计算进程组发信号。
已有输出目录会拒绝覆盖。无样本或采样异常返回 1；参数、平台或附着错误返回 2。
有样本后达到限额返回 0，但必须检查 stop_reason，不能把限额结束当作作业结束。

## 4. 如何用证据决定下一步

目标是回答停顿出现时哪些 rank 在计算、哪些在等待，以及等待是否与存储或内存压力同步。
待检验的假设包括 MCP 复读/回写等待、页回收或缺页、rank 工作量不均和数学库线程争用；
当前没有证据选定其中任何一个为根因。

1. 先验证 rank 覆盖和缺失字段；在可重复输入上比较开关采样的运行，评估扰动。
2. 关联 D/R 线程、wchan/内核栈、read/write 字节、major faults、swap/回收、PSI 和设备计数。
   D 观察次数不是持续时间，更不是“磁盘坏/内存不足”的结论。
3. 若同时需要知道 MATRIX、SETLAG 或 IMPROV 的位置，再实施低频、可关闭的每 rank 内部阶段计时。
   不为计时添加 barrier，不将 collective 的等待时间全部称作传输时间。
4. 证据确定主导阶段后，固定初始轨道、输入与迭代控制，再比较不同 rank 数或针对性修改。
   比较能量、轨道、迭代轨迹和停止状态；普通版/memory 版各自使用自身基线。

本次没有性能测量、原子数值结果或加速结论。预期首轮服务器采样能缩小等待类型和涉及文件范围；
它可能因权限、命名空间或短时事件而不足以定位，届时需补充证据，不能直接改算法归因。

## 5. 本地验证与参考

回归使用临时模拟 `/proc`，覆盖 PID 重用、重新挂靠、读取中退出/替换、权限不足、
rank 冲突、环境白名单、D/R 栈限频、cgroup 祖先限额、线程/fd/输出限额及退出路径。
CTest 注册名为 `rmcdhf_mpi_sampler`。采样器首版通过 15 项 Python 回归。
随后移除普通 MPI 的 NP 入口，在独立 GNU 15 严格检查构建中完成编译与链接；选定 19 项
CTest 通过，包括本采样器。详见[NP 范围变更验证](rmcdhf_nuclear_polarization.md)。未运行科学计算。

```fish
../graspkit-tools/.venv/bin/python test/test_sample_rmcdhf_mpi.py
../graspkit-tools/.venv/bin/python test/npol/check_upstream.py .
```

字段定义：[proc_pid_stat](https://man7.org/linux/man-pages/man5/proc_pid_stat.5.html)、
[内核 PSI](https://www.kernel.org/doc/html/latest/accounting/psi.html)、
[cgroup v2](https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html)、
[磁盘统计](https://www.kernel.org/doc/html/latest/admin-guide/iostats.html)。
