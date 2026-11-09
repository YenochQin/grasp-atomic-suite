# 独立内存 MPI 严格收敛程序

> 状态（2026-11-09 更新）：strict 一致性比较时序已修复；Linux Check/Release 的
> 1/2-rank 重启均在 11 轮满足原严格判据。后续警告修复矩阵通过，全开配置各 93/93；生产体系物理收敛未验证。
> 日期：2026-10-13。复制基准：工作区原版 GRASP 提交
> `488e875d1a7ec631baf4ad33682798bc795b7435` 的 `src/appl/rmcdhf90_mem_mpi`。

最新 [SCNSTY 验证](../reviews/strict_scnsty_normalization_2026-11-09.md)确认：
synthetic Z=10 `1s^2` 原 `6.17e-6` 平台来自未归一化的新解与已归一化旧轨道的幅度差。
仅 strict 将 CONSIS 移至已有归一化之后、阻尼之前；公式、ACCY、能量门槛、连续两轮
要求和测试输入不变。此前 memory IOS 与 100 轮平台证据见
[历史 Linux 报告](../reviews/f2018_linux_davidson_restart_2026-11-09.md)。

新目标是 `rmcdhf_mem_strict_mpi`，源码在
[`src/appl/rmcdhf90_mem_strict_mpi`](../../src/appl/rmcdhf90_mem_strict_mpi)。
它独立于原版 `rmcdhf*` 系列；旧 orbopt 目标已移除。该程序不需要
`GRASP_STRICT_SCF` 环境变量，严格判据始终启用。

复制保留原版内存 MCP、径向求解、阻尼公式、正交化和节点计数算法；
strict 的一致性比较时序修复会改变 SCNSTY 及由它控制的阻尼开关。
适配包括独立目标/程序名、停止策略和输出诊断；网格初始化与参数默认值恢复原版，容量通过 suite 的 parameter_def 转导层提供。
没有引入普通 orbopt MPI 版本的 pair transaction、节点保护或 root tracking。
`libmod`、`lib9290`、`libdvd90`、`mpi90` 使用 suite 的现有共享库，未另建副本。
因此它不是旧服务器二进制的逐位复刻；比较时记录源码版本、编译器、库和参数。

## 停止判据

EOL 每轮结束后同时要求：

1. 所有可变轨道的 `SCNSTY` 有限且非负，最大值 `<= ACCY`。
2. 本轮与前一轮的加权总能量均有限，且
   `abs((E_new - E_old) / E_new) < 0.001 * ACCY`。
3. 上述条件连续两轮成立。任何一次未通过都会将计数归零。

第一次 SCF 迭代不参与能量判据，也不继承重启前的通过次数，因此 EOL 最早
在第三轮结束。零能量或非有限能量不能通过能量条件。
EAL 没有这类对角化能量，只要求连续两轮轨道条件通过。
达到最大轮数时输出 `SCF termination: ITERATION_LIMIT (not converged)`；
只有 `SCF termination: STRICT_CONVERGED` 表示程序的严格条件通过。
最大轮数退出仍沿用上游的正常程序收尾，不能只按进程退出码判定收敛。

实现见 [`strict_scf.f90`](../../src/appl/rmcdhf90_mem_strict_mpi/strict_scf.f90)
及 [`scfmpi.f90`](../../src/appl/rmcdhf90_mem_strict_mpi/scfmpi.f90)。
`ACCY` 仍是共享数值阈值，并非新增的独立能量容差；它还影响积分、插值等过程。
初次诊断建议保持原 ACCY，只改变停止逻辑。例如 `H=0.05` 时，原自适应
`ACCY=H**6=1.5625e-8`，相对能量阈值是 `1.5625e-11`。

`SCNSTY` 保留 `sqrt(UCF) * max(abs(delta P)+abs(delta Q))` 公式；strict 的差值现取自
归一化后的未阻尼建议轨道与旧归一化轨道，发生在阻尼和正交化之前。
它不是阻尼后的保存步长，也不是独立的变分方程残差。2026-11-09 之前的 strict
比较混用了归一化前后的尺度；旧报告中的 SCNSTY 与当前值不能直接按相同语义解读。
低占据轨道仍需检查相邻快照的相位对齐重叠、径向分布、有效节点及目标物理量。
关联轨道的节点不必服从谱学节点数；不要为使节点数相同而将其改设为谱学轨道。

## 输出

`disks` 指定的串行 I/O 目录（无 `disks` 时为启动目录）包含：

| 文件 | 内容 |
| --- | --- |
| `rwfn.out` | 最新一轮径向波函数，保持原 GRASP 格式 |
| `rwfn.out.iter000000` | 读入、插值及正交化后的初始轨道；不是输入文件的逐字节副本 |
| `rwfn.out.iter000001` 等 | 每个外层 SCF 轮次的归一化、正交化后径向快照 |
| `rwfn.out.scf.csv` | 高精度加权能量、相对变化、ACCY、最大 SCNSTY、两项判据及连续通过次数 |
| `rmix.out`、`rmcdhf.sum`、`rmcdhf.log` | 上游格式的混合系数、总结和部分输入记录 |

CSV 的 `T/F` 表示逻辑值；首轮或无法形成能量比较时，`relative_energy_change`
是 `HUGE(double)` 哨兵值，必须结合 `energy_ok=F` 解释。
每轮 CSV 都刷新到磁盘。波函数输出使用 `STATUS='REPLACE'`，避免旧记录尾部残留。
每次续算使用独立输出目录；更短的新作业不会清理旧作业多余的高轮次快照。
快照只包含径向轨道，不能替代完整 `.mix`/CSF 检查点。

## 从最终 .w 续算

1. 建立独立的续算目录，保留原输出。将同一模型的最终 `.w` 复制为串行 I/O
   目录中的 `rwfn.inp`；将对应 AS6 的 `.c`/CSF 文件复制为 `rcsf.inp`，保留
   同一个 `isodata`。不要再次生成 Thomas–Fermi 初猜，也不要把 AS5 波函数覆盖进来。
2. 保持相同的 CSF 内容和顺序、EOL 能级选择及权重。径向文件不会替你恢复这些设置。
   上游入口在启动时重建/对角化 CI 矩阵，不依靠 `rmix.inp` 恢复前一轮本征向量。
3. 保留完整且匹配的 MPI MCP 文件。它们与进程数、rank 及 CSF 顺序有关。
   进程数改变或 MCP 丢失时，应从相同 `rcsf.inp` 重做 `rangular_mpi`，无需重做低层 SCF。
   检查 `disks`：第一行是串行 I/O 路径，后续是各 rank 的临时目录基路径；
   实际 MCP 目录还由共享 MPI 路径代码管理。不能只把文件放在 shell 的当前目录。
4. 显式设置原网格与光速。suite 已恢复原版有限核默认 `N=NNNP=590,NNN1=600`；若旧作业使用
   `N=590,H=0.05`，新输入也应显式指定这些值及原 `RNT,HP,ACCY`。
   读取 `.w` 会插值到本次网格，不会自动采用文件网格。
   编译容量 `NNNP` 与实际 `N` 是两个参数，见[公共参数说明](common_parameters.md)。
5. AS6 外层为 n=9 的例子：保留原更新顺序，只选择第 9 层 13 条轨道可变，
   其余冻结；谱学轨道输入留空。保留 Method 3、原阻尼及更新顺序正交化。
   最大轮数 100 可作为首轮诊断上限；达到上限仍需看残差，不能标记为收敛。
6. 保存实际 stdin/stdout、二进制路径/hash、输入 `.w`/CSF/isodata 的 hash、
   `.sum`、CSV 和各轮快照。检查最后连续几轮轨道形状及节点稳定后，再比较不同模型。

这些是保持模型不变的续算设置，不代表对所有体系都适用的物理精度阈值。
实际 AS6 输入属于计算数据，保留在工作区 `temp`，不加入本仓库。

## 构建及验证

从本仓库根目录，在已初始化编译器、BLAS/LAPACK 和 MPI 的环境下运行：

```sh
export FC=/usr/bin/gfortran-15 OMPI_FC=/usr/bin/gfortran-15
cmake -S . -B build-strict -DGRASP_ENABLE_MPI=ON
cmake --build build-strict --target rmcdhf_mem_strict_mpi test.strict_scf_logic --parallel 4
ctest --test-dir build-strict -R 'strict_scf_logic|grasp_grid_patcher' --output-on-failure
```

上述命令适用于新的构建目录；本次参数恢复已在全新 `build-upstream-parameters` 中全量编译。
参数恢复前的旧构建产物不会自动更新。上述命令的可执行文件在 `build-strict/bin/rmcdhf_mem_strict_mpi`。在 Linux 服务器上重新编译，
不要将 macOS 产物当作 Linux 可执行文件。正常 suite install 也会安装此新目标。

停止策略单元测试覆盖能量平台但轨道未稳、轨道先达标、两轮连续通过、中断归零、
重启计数清零、NaN、零能量、冻结轨道和 EAL；测试调用实际生产模块。
网格补丁器已登记新入口，继续审计集中容量定义、转导接口与各入口初始化。

Linux 端到端测试默认使用同一 suite 构建的准备工具，运行：

```sh
export FC=/usr/bin/gfortran-15 OMPI_FC=/usr/bin/gfortran-15
cmake -S . -B build-strict -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-strict --target rmcdhf_mem_strict_mpi rmcdhf_mem_mpi \
  rwfnestimate rangular_mpi test.strict_scf_logic --parallel 4
ctest --test-dir build-strict -R 'strict_scf_logic|strict_mem_restart_mpi' --output-on-failure
```

注册的 1/2-rank 合成两电子用例检查：轮数上限不会误报收敛、与同库内存基线的
首轮径向更新一致、`.w` 重启读入一致、快照完整，以及最终两轮同时满足判据。
这是小体系工程回归，不替代生产模型的收敛测试。2026-11-09 修复后的 Linux O0/O2
实测单轮、读写和最终两轮判据均通过，1/2 ranks 均于 11 轮收敛，详见文首报告。

此前 macOS 编译成功，停止策略和网格补丁器测试通过。端到端测试遇到继承自
上游的 `systemmem.f90`/`systemfreemem.f90` 依赖 `/proc/meminfo`，无法在 macOS
执行；无该文件时测试明确返回 skip。内存检测算法保持原版，没有用其他求解器
代替并宣称已完成 MPI 验证。
