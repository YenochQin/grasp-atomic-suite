# 内存 MPI 程序的 P0 分阶段诊断

状态：已实现；小型工程回归用于检查记录内容及数值不变性，不代表已修复 Ni 的失稳。
复制基准为原版 GRASP 提交 `fb9a75b320c8d3b67c5c1246fe6d140eb699b856`
中的 `src/appl/rmcdhf90_mem_mpi`。2026-10-13 先复制并逐字节核对 130 个
Fortran 文件和 CMake 文件，再修改 Suite 内的副本。上游源码未修改。

目标 [`rmcdhf_mem_diag_mpi`](../../src/appl/rmcdhf90_mem_diag_mpi/CMakeLists.txt)
独立于 `rmcdhf_orbopt*` 和 `rmcdhf_mem_strict_mpi`。它保留原版的径向求解、
阻尼、正交化、节点约束及停止判据；包括原版“轨道或能量条件通过即可结束”的行为。
本目标没有加入严格收敛、P0 裁剪、自动重试或节点保护。
共享库使用 Suite 的 `libmod`、`lib9290`、`libdvd90`、`mpi90`，不再复制。
因此不能把它称为服务器旧二进制的逐位复刻；比较时需记录编译器、库及运行参数。

## 诊断对象

Ni 松弛测试中，`3d/4d` 的原点展开系数 P0 在能量变化很小、自洽指标仍下降时
增长了多个数量级。现有 stdout 只在一次轨道更新及正交化后输出，不能区分
异常来自解的匹配、归一化、阻尼，还是同对称性轨道正交化。
诊断只读取上述边界上的数据，不作积分、不修正波函数、不增加 MPI 集合操作。

改动所在文件：

- `p0_trace.f90`：CSV 输出与调用上下文；新模块。
- `rscfmpivu.f90`：独立程序名、打开和关闭输出。
- `scfmpi.f90`：初始状态、SCF 轮次和轮末存储轨道。
- `solve.f90`：原始能量估计、内部迭代、齐次/非齐次解匹配及变分修正。
- `improvmpi.f90`：更新入口、SOLVE 返回、归一化、阻尼前后。
- `dampor.f90`：混合之后及阻尼内部再次归一化之后。
- `orthy.f90`：投影对象、参考轨道、投影结果及正交化归一化。

其余复制的 Fortran 文件保持与上述上游提交相同。

## 输出及开关

| `GRASP_P0_TRACE` | 行为 |
| --- | --- |
| 未设置或 `1` | 只输出 rank 0，默认用于诊断 |
| `0` | 关闭诊断，用于同二进制对照 |
| `all` | 每个 rank 分别输出，用于检查并行差异 |

rank 0 在 `disks` 指定的串行 I/O 目录写 `rmcdhf.p0.rank0.csv`。
非零 rank 在 Suite MPI helper 提供的并行工作目录写各自的
`rmcdhf.p0.rank<rank>.csv`；`permdir` 只在 rank 0 有效，未修改共享 MPI helper。
每条记录写入后刷新。正常结束会关闭文件；原版 `STOP` 前已刷新记录仍可读取。
输出失败会打印 `[P0_TRACE]` 警告并关闭本 rank 的诊断，不改变求解路径。
每次运行替换同名诊断文件，因此应使用独立作业目录。

CSV 中包含：

- `cycle/update/inner_iteration`：SCF 轮次、本轮轨道更新序号、SOLVE 内部轮次。
- `active_orbital`：触发这次更新的轨道。`orbital/label/n/kappa` 标识实际记录的轨道；
  ORTHY 可能同时改变另一个同 kappa 轨道，二者不能混为一谈。
- `peer`：正交投影的另一个轨道编号，无另一个轨道时为 0。
- `p0/energy/gamma/scnsty/method/mtp`：当前系数及求解上下文。
- `p2..p7/q2..q7/r2..r7`：前六个非零网格位置上的原始数组值及半径。
  超出本次有效表长的 P/Q 样本填 0。P0 是展开系数，并非物理上的 P(0)。
- `overlap/scale/norm/mixing/alpha/match_denominator`：下表所述阶段的辅助量。
  不适用于某条记录的辅助列填 0，不能将其解释为实际测得的零范数等。

| stage | 数据来源及辅助量 |
| --- | --- |
| `initial`、`cycle_end` | 全部存储轨道；`initial` 是加载/插值/正交化后的初值 |
| `update_begin`、`solve_iteration_begin` | 存储的 PZ/PF/QF；当前更新及内部迭代入口 |
| `energy_estimate` | `EIGEN(J)` 原始估计，发生 EPSMIN 截断之前；P0/PQ 仍为存储值 |
| `before_match` | 非齐次试探解；`match_denominator=QJPIH-QJPOH` |
| `matched_solution` | 匹配后的试探解；记录 `alpha` 与匹配分母 |
| `variational_solution` | Method 2/4 的修正解；`scale=DELE`，`alpha` 为变分匹配系数 |
| `solve_return` | 成功返回的 P0/P/Q；尚未执行 IMPROV 的归一化 |
| `solve_failure` | 失败时的存储轨道；不能误认为失败试探解 |
| `normalized_solution` | IMPROV 归一化后的试探解；`norm` 为归一化前值，`scale=DNFAC` |
| `before_damp`、`after_damp` | 存储轨道的阻尼前后状态；`mixing=ODAMPJ` |
| `damp_mixed`、`damp_normalized` | 混合后的存储状态及内部归一化；记录 `mixing/norm/scale` |
| `orth_begin`、`orth_end` | 每个被正交化轨道的入口和尾部裁定后的出口 |
| `orth_before_projection`、`orth_reference`、`orth_after_projection` | 投影的目标、参考及结果；`peer/overlap` 标识配对 |
| `orth_normalized` | 投影后的归一化；`norm` 为归一化前值，`scale` 包含程序的符号调整 |

建议按 `(cycle, update)` 追踪 `3d/4d`，从 `matched_solution` 到
`normalized_solution`、`damp_mixed`、`damp_normalized`、`orth_*` 查找首次异常增长。
必须同时看第一批径向值、匹配分母和归一化因子，不能单独用一个 P0 阈值判定物理错误。

## 构建与现有作业脚本

在服务器的 Suite 仓库中构建：

```sh
cmake -S . -B build-p0-diag -DGRASP_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-p0-diag --target rmcdhf_mem_diag_mpi --parallel 4
```

产物在 `build-p0-diag/bin/rmcdhf_mem_diag_mpi`。无需覆盖原版可执行文件。
Debug 用于工程诊断；与既有生产结果对照时，应同时记录并尽量匹配原编译选项，
不要假设不同编译器或优化级别会在完全相同的轮次失稳。
共享参数仍由 `suite_parameters_M.f90` 提供；本次未改变 NNNP、N、H 或 ACCY。
新程序的网格入口已登记到 `atomic-suite` 网格补丁布局。

保持原 CSF、核参数、初始波函数、目标能级、更新顺序、阻尼和 MPI 进程数，
通过现有作业脚本的 `GRASP_RMCDHF_MPI` 覆盖选择新程序：

```sh
export GRASP_RMCDHF_MPI=/path/to/grasp-atomic-suite/build-p0-diag/bin/rmcdhf_mem_diag_mpi
export GRASP_P0_TRACE=1
# 然后按原作业脚本的参数提交到独立目录。
```

这不是提交授权。对于当前保存结果的脚本，原生诊断 CSV 不匹配带 AS 前缀的归档规则，
通常会保留在作业 `.work/rmcdhf.p0.rank0.csv`，请和求解器 stdout 一并取回。
失败的 Ni checkpoint 不应覆盖正常 baseline；诊断重放应从原先验证过的相同输入开始。

## 验证

```sh
cmake -S . -B build-p0-diag -DGRASP_ENABLE_MPI=ON \
  -DGRASP_PREPARATION_BINDIR=/path/to/grasp/bin
cmake --build build-p0-diag --target rmcdhf_mem_diag_mpi \
  test.p0_trace_projection test.shared_grid_capacity --parallel 4
ctest --test-dir build-p0-diag \
  -R 'p0_trace|shared_grid_capacity|grasp_grid_patcher' --output-on-failure
```

`p0_trace_projection` 调用实际 ORTHY 和诊断模块，用可精确计算的二维内积
检查 P0 投影、归一化、参考轨道记录、rank 开关、NaN 保存及输出失败处理。
`p0_trace_mpi_1/2` 使用合成两电子模型，对照关闭、仅 rank 0 和全部 rank 诊断，
检查各阶段记录以及 `rwfn.out/rmix.out` 的逐字节一致性。
MPI 测试需要 Linux `/proc/meminfo`、准备程序及可用的本地 MPI 通信环境。
这些测试不是 Ni 生产失稳的复现，也不证明原求解器的物理正确性。

2026-10-13 本地 Debug 构建通过，五项针对性 CTest 通过：共享容量、网格补丁器、
投影追踪及 1/2-rank MPI 对照。MPI 需要允许本地 socket 通信；首次沙箱运行因
该限制失败，在允许本地 MPI 通信的环境完成验证。另有原版 Fortran 接口兼容性
警告，未改写这些继承代码。未提交或运行 Ni 生产作业。
