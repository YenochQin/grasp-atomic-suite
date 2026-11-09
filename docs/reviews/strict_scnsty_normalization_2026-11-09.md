# Strict 重启 SCNSTY 平台：归一化时序修复

后续更新：[32 配置矩阵与警告基线](f2018_matrix_warning_baseline_2026-11-09.md)已通过；
新增设施测试后全开配置各 75/75。本文保留 strict 修复阶段的 74/74 和数值证据。

日期：2026-11-09。基线为 `f395146` 加
[上一轮 Linux 修复](f2018_linux_davidson_restart_2026-11-09.md)的未提交工作树。
本报告保留其 71/73 和旧版 100 轮未收敛记录，不覆盖旧失败证据。

**已定位并修复 strict 的比较时序；原 1/2-rank 重启用例均在 11 轮满足严格判据。**
没有改 ACCY、能量阈值、连续通过次数、输入、阻尼参数或最大迭代次数。
SCNSTY 的比较对象由未归一化的新解改为已归一化、尚未阻尼的新建议轨道；
这是有意的语义修正，不能描述成纯 F2018 语法变换。

## 实验与根因证据

目的：区分真实轨道形状停滞、归一化幅度差，以及重启/尾部数据问题。
先用原 `strict_mem_restart_mpi_1` 复现 100 轮后收敛断言失败，再在 strict IMPROV
现有归一化前后记录同一更新的范数、SCNSTY、归一化形状差及纯幅度差预测。

输入保持 finite-nucleus synthetic Z=10 `1s^2`，
`N=NNNP=590, NNNW=127, RNT=2e-7, H=0.05, HP=0, ACCY=1.5625e-8`；
Method 3、单目标 EOL，从原首次一轮输出重启，最大 100 轮。
GNU 15.3.0、单 rank，BLAS/OpenMP 均单线程。若是幅度问题，归一化后的形状差
应远小于原 SCNSTY；若是实际形状误差，该差仍应停在相同量级。

诊断末轮得到：

| 量 | 修改前第 100 轮 |
| --- | --- |
| 未归一化建议轨道的范数平方 | `0.99999630588248001` |
| 原 SCNSTY | `6.1662534234748093e-6` |
| 归一化建议轨道与旧轨道的加权 P/Q 最大差 | `7.9167005690195972e-14` |
| 纯幅度差预测 | `6.1662534779606720e-6` |

纯幅度差预测使用
`sqrt(UCF)*abs(sqrt(DNORM)-1)*max(abs(P_old)+abs(Q_old))`，
与原 SCNSTY 相差约 `5.45e-14`。因此该平台由比较双方尺度不一致主导。
只加观测的诊断版与原失败运行的全部轨道及 CSV 逐字节一致。

调用顺序原为 `SOLVE → QUAD → CONSIS → normalize → damping`。
`PF/QF` 是前一轮已经归一化的轨道；Method 3 的径向建议解尚待调用者归一化。
`CONSIS` 提前读取 P/Q，把归一化前的幅度偏差混入严格轨道条件。形状已稳定时，
这个偏差仍阻止连续两轮通过。原版较早的能量停止可能使这一问题不表现为测试失败。
本轮并未把原始 DNORM 的离散误差来源进一步分解为积分、估能或边界误差。

## 修复范围与不变条件

只在[严格 IMPROV](../../src/appl/rmcdhf90_mem_strict_mpi/improvmpi.f90)中，
将已有 `CALL CONSIS(J)` 移至已有 P0/P/Q 归一化之后、DAMPOR 之前。
`CONSIS` 的 `sqrt(UCF)*max(abs(delta P)+abs(delta Q))` 公式和原网格区间不改。
相应更新 CONSIS/STRICT_SCF 注释，明确它测量归一化的**未阻尼建议轨道**变化。

- 仍要求所有可变轨道 `SCNSTY <= ACCY`，EOL 能量相对变化 `<0.001*ACCY`，连续两轮成立。
- 首轮没有有效能量比较；重启从 streak=0 开始。`ITERATION_LIMIT` 不计作收敛。
- 比较在阻尼前完成；强阻尼导致保存步长变小不能触发虚假通过。
- 未改变 SOLVE、积分、核模型、根选择、网格、初始轨道、Method 或输入阻尼参数。
- 普通、memory、diag、CSFG 求解器的 CONSIS 调用时序保持原状；不在本次修复中统一其停止策略。
- strict 的输出 SCNSTY 和由它控制的阻尼开关可能改变；这是修复的行为变化，不要求复现旧假平台。

这使两次被比较的轨道处于同一归一化尺度，不代表原始径向方程残差已小于 ACCY。
SCNSTY 仍是轨道更新指标，不能将本回归视为连续方程或物理模型精度认证。

## 回归设计

[strict_orbital_consistency.f90](../../test/f2018/strict_orbital_consistency.f90)
编译真实 strict IMPROV、CONSIS、DAMPOR，并使用真实 QUAD/RINT。
仅固定径向方程装配与 SOLVE 输出，测试的是生产调用顺序，不是重写一份归一化算法。
解析夹具为 `[0,9]` 上 `P=r/sqrt(243), Q=0`，范数积分恰为 1；
完整 Newton–Cotes 面板对该二次被积函数精确。

四个案例覆盖：2 倍幅度、1/4 倍幅度、真实 P/Q 形状旋转、原 ACCY 以下的小形状变化。
用旋转角的解析表达式验证 SCNSTY，真实变化案例另施加强阻尼，使保存后的步长小于
ACCY，仍要求未阻尼 SCNSTY 大于 ACCY。既检查缩放不产生假不收敛，也检查阻尼不产生假收敛。
修改前在第一个纯缩放案例报 `Normalized proposal consistency`；修改后通过。

新增 `f2018_strict_orbital_consistency` 已加入 required 清单。
原 `strict_scf_logic` 和 `strict_mem_restart.py` 均未改判据/断言；后者继续验证首轮等价、
重启 I/O、快照完整性、有限值、最后两轮轨道及能量条件。

## 收敛结果与验收

同一输入在 Check/O0、Release/O2 的 MPI 1/2 ranks 均于 11 轮重启迭代后
输出 `STRICT_CONVERGED`。归档的单-rank O0/O2 最后两轮完全相同：

| 轮次 | max SCNSTY | 相对能量变化 | pass streak |
| --- | --- | --- | --- |
| 10 | `2.4537930189302858e-9` | `5.0079772408038271e-13` | 1 |
| 11 | `6.4473567455537642e-10` | `2.3225401696476121e-13` | 2 |

两列分别低于原门槛 `1.5625e-8` 和 `1.5625e-11`。一轮限定运行仍报告未收敛，
与普通 memory 的首次径向更新对照仍通过。归档 O0/O2 的 12 份快照、最终波函数、
CSV 和 mixing 共 **15 文件逐字节一致**，失败/修复两侧的输入 hash 一致。

相对旧版 100 轮平台，修复后最终加权能量变化为 `+1.91420213013771e-11 Ha`，
P/Q 最大逐点差分别为 `5.197537183573111e-10`、`5.760736287498161e-11`；
轨道能量变化 `-2.3590018827235326e-12 Ha`。这些是实测软件对照，旧平台不是物理 oracle。

Linux GNU 15.3.0 全量结果如下，serial/MPI/CSFG/OpenMP 均开启：

| 配置 | 实际编译命令 | 全量构建 | CTest | 强制验收 |
| --- | --- | --- | --- | --- |
| Check/O0 | 3,254 | 通过 | 74/74，零失败/跳过 | 返回 0，`strict_acceptance=true` |
| Release/O2 | 3,254 | 通过 | 74/74，零失败/跳过 | 返回 0，`strict_acceptance=true` |

包含原停止逻辑、两项重启测试、新生产调用顺序回归、NP 15/15、生产 CI 阻尼 14/14、
Davidson 8/8、jj2lsj 3/3 和全部当前已注册消费者。原版来源保护仍为 492/492。
没有扩大编译放行或改变 fail-closed 验收规则。
**这是 required.json 当前范围通过；完整 P0–P9、配置矩阵和 F00–F10 仍未完成。**

## 复现与证据位置

从 suite 根目录运行；工具链沿用上一轮独立 Linux 构建目录：

```sh
export FC=/usr/bin/gfortran-15 OMPI_FC=/usr/bin/gfortran-15
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-linux-check --jobs 4
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-linux-release --jobs 4
ctest --test-dir build-gf15-linux-check \
  -R '^(f2018_strict_orbital_consistency|strict_scf_logic|strict_mem_restart_mpi_[12])$' \
  --output-on-failure -V
```

- `build-gf15-linux-check/audit/strict-platform/`：修改前源文件、旧 Check/Release 验收摘要、
  `unit-before.log`/`unit-after.log`、`restart-after.log`、`fixed-comparison.json`。
- 同目录的 `diagnostic-run/`：只加观测的 100 轮结果，`diagnostic.log` 保留幅度/形状分解；
  临时 `[DEBUG-scnsty]` 输出已从源码移除。
- `build-gf15-linux-{check,release}/audit/strict-platform-fixed/`：修复后输入、完整运行日志、11 轮快照及 CSV。
- `build-gf15-linux-{check,release}/audit/acceptance/`：最终注册清单、编译命令检查、全量构建、CTest XML 及源码 hash。
- 上一轮 `audit/strict-restart-unconverged/` 和 baseline 的归档仍保留原失败证据。

没有 install、提交、push 或 SLURM 作业。正常数值停止已经在此软件夹具上验证；
真实多根、实际原子属性/跃迁和物理收敛仍按迁移计划单独验收。
