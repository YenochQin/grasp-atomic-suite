# 公共参数集中管理

> 文档状态：已实施、现行配置说明。
> 核查日期：2026-10-05；源码基准：`8b298a6`。
> 公共参数和共享网格初始化已接入六个入口；工程测试与物理收敛是不同证据。当前使用入口见本文及[网格实验方案](../plans/grasp_grid_parameter_tuning.md)。
> 全部文档状态见[分类索引](../README.md)。

所有串行、MPI、内存版本及 RHFS/RIS/g_J 使用同一套公共默认值。
修改入口只有 [`suite_parameters_M.f90`](../../src/lib/libmod/suite_parameters_M.f90)。
运行时变量继续保存在 `grid_C`、`def_C`、`coun_C`，由共享模块
`radial_grid_defaults` 初始化。`parameter_def` 只转导旧参数名，没有第二份定义。

## 可调整参数

| 参数 | 当前默认值 | 含义 |
| --- | --- | --- |
| `NNNP` | 2990 | 径向数组编译容量；旧求解路径也用于尾部边界与连接点回退 |
| `NNNW` | 127 | 轨道数组编译容量 |
| `KEYORB` | 215 | 原有整数打包参数，需保持打包范围约束 |
| `FINITE_N` | `NNNP` | 全部有限核程序的默认实际点数 |
| `FINITE_RNT_SCALE` | `2D-6` | 有限核 `RNT=FINITE_RNT_SCALE/Z` |
| `FINITE_H` | `0.05D0` | 有限核默认指数步长 |
| `POINT_N` | `MIN(220,NNNP)` | 点核默认实际点数 |
| `POINT_RNT_SCALE` | `EXP(-65D0/16D0)` | 点核 `RNT=POINT_RNT_SCALE/Z` |
| `POINT_H` | `0.0625D0` | 点核默认步长 |
| `DEFAULT_HP` | `0D0` | 两种核模型的公共网格参数 |
| `DEFAULT_ACCY` | `0D0` | 0 表示自适应 `H**6`；正值表示固定默认阈值 |
| `NODE_THRESHOLD` | `0.05D0` | 四个轨道程序的默认节点振荡筛选阈值 THRESH |

`NNN1=NNNP+10`、`NNNWM1=NNNW-1`、`NNNWM2=NNNW-2` 自动派生，不能分别修改。
原来性质程序的有限核 590 点默认值已按用户要求统一为 `NNNP`。
数值默认值集中后仍需通过物理量的网格收敛判断合适的 H/N/RNT；统一配置本身
不等于已经达到计算精度目标。

当前 `N=2990,H=0.05` 不宜直接作为通用网格起点：审查已在原版
`DCBSRW` 最小径向诊断中复现非有限值，不能以远端最终清零排除中间溢出。
本次文档修订保留既定 Fortran 默认值；在独立测试副本上先显式选择并应用
规程中的试验配置，例如 D2 的 `NNNP=2990,N=1965,H=0.015,ACCY=1e-10`，
再检查各阶段有限性与物理收敛。该例是候选，不是已认证生产参数；
诊断范围与复现过程见[深度审查报告](../reviews/grasp_grid_methodology_review.md)。

原值与候选值、控制变量试验、测试数据和误差预算见
[径向网格参数调整与验证规程](../plans/grasp_grid_parameter_tuning.md)。

物理常数继续由现有 `SETCON` 统一设置；特定轨道算法的实验开关继续由
`orbopt_control` 管理，不混入公共网格参数。

## 初始化与覆盖顺序

1. 读取核参数，获得 NPARM 和 Z。
2. 六个入口统一调用 `SET_RADIAL_DEFAULTS(NPARM,Z)`，设置 RNT/H/HP/N/ACCY。
   g_J 与 RHFS 复用同一入口；这六处覆盖全部 10 个程序。
3. 保留程序已有的交互输入，允许修改网格或光速。
4. 用最终 H 调用 `UPDATE_RADIAL_ACCURACY`。默认 ACCY 为 `H**6`，固定配置时
   保持 `DEFAULT_ACCY`。RMCDHF 在提示用户是否修改 ACCY 前执行此步骤。
5. 用户随后显式输入 ACCY。MPI 版本保留原有广播顺序。
6. `VALIDATE_RADIAL_GRID` 检查容量、有限正数和 HP 非负性，然后生成网格。
   验证不会覆盖用户已输入的 ACCY。
7. 轨道文件读取通过 `INTRPQ` 插值到当前计算网格，不将初始化参数替换成
   文件中轨道的原网格。`rmcdhf_orbopt_mpi` 原有的 `GRASP_COUNT_ACCY` 和
   `GRASP_COUNT_THRESH` 环境覆盖仍在后续原位置应用。

`GRASP_COUNT_ACCY` 实际覆盖全局 ACCY，不只影响节点计数；前面的轨道
初始化/插值可能已经使用不同阈值。网格试验应清除覆盖，或记录各阶段取值。
`GRASP_STRICT_SCF` 只按所用变体的实际实现解释；原版的收敛退出逻辑与
suite 普通 MPI 严格模式不同。统一按规程检查变分轨道残差、目标量和有效尾部。
容量审计、全量编译、实际执行路径和物理验收是不同层次的证据。

## 修改方法

直接修改中央 Fortran 文件即可。也可以使用原有脚本预览：

```sh
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --layout atomic-suite --grasp . --nnnp 2990 --n 1965 --h 0.015 \
  --rnt-scale 2e-6 --hp 0 --accy 1e-10 --diff
```

未显式指定 `--n` 时，`FINITE_N=NNNP` 会自动跟随容量。显式 `--n` 会为所有
有限核程序设置相同的独立实际点数，必须不超过容量。
`--point-n/--point-h/--point-rnt-scale` 调整点核默认值。
`--accy 1e-10` 设置固定阈值，`--accy 0` 恢复自适应阈值。
示例是接口用法，不是对某个元素的推荐数值。

追加 `--apply` 才会写入；脚本只修改中央文件，并审计全源码中重复的容量定义、
脱离公共初始化的径向默认赋值和缺失的入口。应用/检查/备份/恢复的具体说明见
[scripts/README.md](../../scripts/README.md)。完整原版 GRASP 的默认布局保持原有补丁行为。

修改后在已开启所需串行/MPI 目标的构建目录全量编译、测试并安装：

```sh
cmake --build build-all --clean-first --parallel 4
ctest --test-dir build-all --output-on-failure
cmake --install build-all
```

确认调用路径指向本次安装的程序，并查看输出记录中的实际网格参数。
源码树中的旧 `.mod` 会导致混合容量；CMake 的检查及共享数组容量测试继续启用。

## 集中参数改造时的工程验证（历史记录）

- 10 个应用成功编译；串行与 MPI 使用同一配置模块。
- MPI 构建通过 11 项 CTest，纯串行构建通过 10 项。
- 网格脚本通过 16 个离线用例，包括重复声明、局部同名变量、固定/自适应
  ACCY 切换、容量缩小时实际点数冲突以及恢复。
- 真实源码副本中仅修改中央文件，随后编译验证 `NNNP=1990`、`NNN1=2000`、
  数组长度、有限核 N/H、点核 N/H/RNT、固定 ACCY 和显式用户覆盖全部生效；
  应用/检查/恢复通过。
- `data/test` 的串行和双进程 MPI g_J 输出与调整前逐字节一致。
  该结论只覆盖这组输入，不代表所有体系增加 N 都不影响结果。
- 大型 RMCDHF 原子计算和集群作业仍需使用外部输入另行做物理验证。
