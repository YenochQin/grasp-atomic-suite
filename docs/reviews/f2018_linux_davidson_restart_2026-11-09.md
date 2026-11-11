# F2018 Linux 验收、Davidson 早收敛与 memory 路径修复

> 2026-11-11 状态更新：文中 strict/diag 独立求解器已按用户决定[移除](../implemented/rejected_rmcdhf_variants.md)；相关目标、命令和测试数量仅作对应版本的历史记录。

后续更新：本文保留本轮 71/73 与原始平台证据。随后已完成 strict 归一化时序修复，
Check/Release 当前清单各 74/74 通过，见[后续报告（已移除）](../implemented/rejected_rmcdhf_variants.md)。

日期：2026-11-09。源码基线：`f395146`（包含严格接口修复 `b8d6faf`）；
结果对应本轮工作树。接续[此前 macOS 验证](f2018_strict_interfaces_npol_2026-11-09.md)，
更新[迁移计划](../plans/f2018_migration_gfortran15_2026-11-08.md)。

**Linux GNU 15.3.0 Check/Release 全量构建通过；各 73 项测试为 71 通过、2 失败、
0 跳过。两个 strict 重启测试仍未通过，强制验收返回 1，`strict_acceptance=false`。**
本轮未完成 P3/P5/P9，也未证明原子物理精度或活性空间收敛。

## 核查与工具链

开始时工作树干净。计划记录的严格接口、NP 路径、SOLVE 意图等修复已包含在
`b8d6faf`；此前报告中“未提交工作树”是当时状态。已知 MPI Davidson 早收敛分支
仍可复现；四项此前因缺少 Linux `/proc/meminfo` 而跳过的测试在本机实际运行。

环境为 Linux x86_64、GNU Fortran 15.3.0、CMake 4.4.5、Open MPI 5.0.11、
OpenBLAS 0.3.34 LP64。`FC/OMPI_FC` 均显式选择 `/usr/bin/gfortran-15`，
MPI wrapper 后端亦为该编译器。非线程专项固定 `OMP_NUM_THREADS=1`、
`OPENBLAS_NUM_THREADS=1`。系统 MPI 的运行探针通过不等于已取得其完整绑定兼容证明。

Check/O0 和 Release/O2 使用独立新目录，serial/MPI/CSFG/OpenMP 全开，
各检查 **3,241 条编译命令**。严格标准、接口、浮点选项和原数值门槛未放宽。
没有全局静态存储、默认 kind 扩宽、参数放行或额外计算线程。
编译器、数学库、参考可执行程序的 hash、实际命令与源码工作树 hash 留在忽略目录中。

## 已实施的三项缺陷修复

### CSFG MPI Davidson 的早残差收敛（F 类）

原常规三对角测试通过，不能覆盖独立 2×2 不变子空间。
新增输入为 6×6 实对称矩阵：左上块 `[[2,1],[1,2]]`，其余对角 `4,5,6,7`，
从第一坐标向量出发，`LIM=6`（保留 GDVD 的重启工作区要求）。解析最低本征值为 1。
同一输入串行通过，修改前 MPI 单 rank 在 MGS_NRM 处发生 SIGFPE。

调试器记录：`NLOOPS=1, KPASS=2` 时直接残差平方 `SQRES=0`，
`CRITR=1e-12`，随后进入归一化，`dnm=0`。`NEWVEC` 的 `NLOOPS>2` 条件
拒绝已经满足残差门槛的解，导致继续预条件化并除零。此结论来自实际断点，
不是仅凭上一轮怀疑或最终能量推断。

只移除[NEWVEC](../../src/lib/libdvd90mpi/dvdson.f90)的直接残差分支轮数限制。
残差使用原 `CRITR`；`TSTSEL` 对能量差及系数的三轮限制、预条件化、归约和重启
策略保留。这是有解析 oracle 的终止分支缺陷修复，早收敛输入的行为有意改变，
不要求复现旧崩溃；不能描述成纯语法等价改写。

局部警告登记：`src/lib/libdvd90mpi/dvdson.f90` 的 NEWVEC 为保持已有接口仍接收
`NLOOPS`，本轮修复后该形参不再读取，新增 `-Wunused-dummy-argument` 告警。
未抑制它；后续接口现代化时可连同调用者/`newvec_I` 一起移除。
本轮没有完成全仓库逐文件警告基线，其余历史告警仍需独立审计。

保留原 6×6 三对角案例，新增串行及 MPI 1/2/4-rank 早收敛测试，全部加入
[required.json](../../test/f2018/required.json)。两种构建中常规/早收敛各 4/4 通过；
检查有限值、本征值、显式残差、归一化，门槛均为 `1e-12`，设置超时。
覆盖生产 `DDOTMPI/DGEMVMPI/NEWVEC/MGS_NRM`，不等于真实多根 CSFG 原子链验收。

### 普通 RANGULAR_MPI 空 rank（F 类）

原 1-CSF/2-rank 准备路径在 `SETSDA` 访问未初始化的 `IC`，严格运行时捕获
IENDC 越界。空 rank 的 `NNONZ=0`，读入循环不会设置 `IC`。
改用同过程已初始化为零的 `ICLAST`；非空路径它等于末列 IC，文件协议不变。

[独立回归](../../test/f2018/sparse_empty_rank.f90)调用真实 SETSDA，依次检查有工作
和无工作的 rank，解码 MCP header、计数、列尾/行号及空记录。
仅该测试副本使用 `-finit-integer=-999999 -fcheck=all`，使旧错误稳定表现为负索引越界；
该诊断选项没有施加于生产目标。修改前失败、修改后 Check/Release 通过。
实际 2-rank P0/strict 路径亦已越过准备阶段。

### Memory SETHAM 的无效 IOS 判断（F 类）

普通 memory、memory MPI、strict、diag 四个 SETHAM 保留了旧文件读取后的
`IOS` 判断，但已改为调用 `read1_mem` 从缓存数组取标签；该接口没有文件 I/O
或 IOSTAT 输出，局部 IOS 从未赋值。GNU 15 的未初始化警告与 Linux 单 rank
`STOP IOS .NE. 0 when reading LAB, NCONTR` 复现相互印证。

移除这四份实现的 IOS 声明与对应判断；保留实际磁盘加载阶段的 I/O 检查、
缓存索引更新、积分计算及矩阵装配顺序。没有用 `IOS=0` 或全局默认初始化掩盖问题。
原版两个受保护文件的可逆差异追加到
[兼容性清单](../reference/rmcdhf_f2018_compatibility.json)，原始来源 hash 不变，
`check_upstream.py` 仍为 **492/492**。

## 完整回归与未通过项

| 检查 | Check/O0 | Release/O2 |
| --- | --- | --- |
| 全量构建、清单及实际编译命令 | 通过，3,241 条 | 通过，3,241 条 |
| 完整 CTest | 71 通过 / 2 失败 / 0 跳过 | 71 通过 / 2 失败 / 0 跳过 |
| 常规 + 早收敛 Davidson | 8/8 | 8/8 |
| 空 rank 独立回归 | 1/1 | 1/1 |
| P0 trace MPI 1/2 ranks | 2/2 | 2/2 |
| NP（含来源保护） | 15/15 | 15/15 |
| 原生产 DDOT/CI 阻尼分支 | 14/14 | 14/14 |
| strict restart MPI 1/2 ranks | 0/2 | 0/2 |

修改前 Linux 全量 Check 为 64/68 通过、4 失败（不是 macOS 的 4 跳过）。
修复后 P0 tracing 的 off/root/all 输出不变性通过；strict 测试中单轮上限、
首次 CSV、波函数重启读写及与原版 memory 第一轮对照通过，但后续 100 轮
**未达到 STRICT_CONVERGED**。完整测试仍记失败，未删除断言、加大轮数或放宽 ACCY。

## Strict 停滞的独立旧新对照

目的：判断新暴露的 strict 未收敛是本轮迁移回归，还是基线已有行为。
如果本轮只修复无定义状态读取，预期相同已定义路径保持一致；不能预先断言会收敛。

输入为现有有限核 synthetic Z=10 `1s^2`，`N=590, NNNP=590, NNNW=127`，
`RNT=2e-7, H=0.05, HP=0, ACCY=1.5625e-8`，从首次一轮波函数重启，最多 100 轮。
参考为 `f395146` 的独立 Git 归档、GNU 15 legacy/O0（保留 `-fno-automatic`），
候选为本轮严格 Check/O0；同数学库、单线程、MPI 单 rank。

实际：两者均运行满 100 轮，终止于 `ITERATION_LIMIT`。全部 **104 个输出文件逐字节一致**：
SCF CSV、101 份逐轮波函数、最终波函数及 mixing。末轮共同为：

| 量 | 旧基线与候选 |
| --- | --- |
| 加权能量 / Ha | `-93.982757016713137` |
| 相对能量变化 | `3.0241408458947123e-16` |
| max SCNSTY | `6.1662534234748093e-6` |
| ACCY | `1.5624999999999999e-8` |
| orbital_ok / energy_ok / pass_streak | `F / T / 0` |

因此该案例的停滞不是本轮修复引入。能量条件成立不能替代轨道条件；
`SCNSTY` 的来源及为何形成此平台尚未定位，本轮不改变其定义或 strict 判据。
两项失败继续阻止验收。旧基线多轮对照只做了单 rank/O0；候选完整测试覆盖 1/2 ranks、O0/O2。

## 普通与 memory 串行一轮对照

目的为确认本轮修复和严格存储设置下的固定一轮软件行为，预期同优化级别旧新完全一致。
输入沿用上述有限核单组态模型、单线程，最大一轮，不作正常收敛声明。
`f395146` 归档分别以 legacy/O0 和 legacy/O2 构建；两者均显式关闭 fast-math 与 FMA contraction。
候选分别采用 Check/O0 与 Release/O2，未用 O0 参考替代 O2。

默认完整两个 solver（`rmcdhf`、`rmcdhf_mem`）在每种优化级别各重复两次，
共八组 solver 对照。每组七个二进制文件及全部科学文本/诊断均一致，
`full_serial_scope=true`，解码与完成状态检查通过，未扩大比较器过滤规则。
同配置重复运行的二进制、科学文件和计算数值稳定；memory 的 stdout 在两次运行之间
仅 Linux 空闲内存报告变化（此次显示为 17.27→17.41 Gb），差异原文保留。
它不支持“跨任意时刻整个 stdout 恒定”的声明，也不能用来忽略分配失败或其他诊断。

## 复现与证据

从 suite 根目录执行，使用工作区共享 Python；MPI 测试需要本机进程通信权限。
两个目录均按计划配置 GNU 15、serial/MPI/CSFG/OpenMP 全开：

```sh
export FC=/usr/bin/gfortran-15 OMPI_FC=/usr/bin/gfortran-15
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-linux-check --jobs 4
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-linux-release --jobs 4
# 上述两次仍预期返回 1：strict restart 的收敛断言失败。
ctest --test-dir build-gf15-linux-check \
  -R '^(f2018_davidson|f2018_sparse_empty_rank|p0_trace_mpi)' --output-on-failure
../graspkit-tools/.venv/bin/python test/f2018/compare_first_sweep.py \
  --reference-bindir build-gf15-linux-baseline/check/bin \
  --candidate-bindir build-gf15-linux-check/bin \
  --output build-gf15-linux-check/audit/first-sweep-new
```

本地忽略目录中的证据：

- `build-gf15-linux-{check,release}/audit/acceptance/`：实际命令、构建、CTest XML、源码 hash 与失败摘要。
- `build-gf15-linux-check/audit/before-linux-fixes/`：初始 Linux 64/68 及四项失败。
- `build-gf15-linux-check/davidson-early-before.log`、`davidson-before-gdb.log`：SIGFPE 与断点证据。
- `build-gf15-linux-check/sparse-empty-before.log`：空 rank 的确定性越界反例。
- `build-gf15-linux-{check,release}/audit/first-sweep-{1,2}/`：两种优化级别、两个 solver、重复两次的输入/输出/hash。
- `build-gf15-linux-baseline/`：归档来源及 O0/O2/单-rank legacy 构建；`audit/toolchain.json` 记录工具链和参考二进制身份。
- `build-gf15-linux-baseline/audit/first-sweep-repeat.json`：重复运行间的内存报告差异。
- `build-gf15-linux-check/audit/strict-restart-unconverged/` 和
  `build-gf15-linux-baseline/audit/strict-restart-1/`：未收敛的候选/基线输入、日志和全部逐轮输出。
- `build-gf15-linux-baseline/audit/strict-restart-comparison.json`：末轮数值及 104 个精确一致文件的 hash。

下一步优先独立定位 strict 的 SCNSTY 平台，保留原失败样例与判据；并继续手写接口/状态
审计、警告基线、统一 case manifest、配置矩阵、真实多根/RCI/跃迁及 `mpi_f08`/大计数。
本轮没有 install、提交、push 或 SLURM 作业，没有恢复任何退役求解器。
