# 参数转导层与原版默认值

> 文档状态：现行实现；2026-11-16 共享容量使用语义名称，原版默认值和网格初始化保持一致。
> 参考源码：工作区 GRASP `488e875d1a7ec631baf4ad33682798bc795b7435`。
> 移除 CSFG 后，网格脚本覆盖十五个普通网格入口与两个 RCI 续算 ACCY 入口。
> 容量由集中配置定义，经 grasp_parameters 转导；历史实验数值不代表当前默认值。

## 参数归属

[`compile_time_capacities.f90`](../../src/lib/grasp_shared/compile_time_capacities.f90) 是容量常量的唯一来源；
[`grasp_parameters.f90`](../../src/lib/grasp_shared/grasp_parameters.f90) 用 `USE compile_time_capacities, ONLY: ...`
直接转导六个语义名称。共享库和各串行/MPI 程序继续 `USE grasp_parameters`，
修改配置不需要逐个修改调用方。配置模块默认值与参考原版相同。
网格的 `radial_grid_defaults` 集中初始化模块不使用，各入口保留原版初始化。

| 配置模块和应用使用的名称 | 默认值 | 含义 |
| --- | --- | --- |
| `orbital_index_radix` | 215 | 轨道整数打包参数 |
| `max_radial_points` | 590 | 径向数组编译容量；部分旧径向路径也将它用于尾部边界 |
| `radial_workspace_capacity` | 600 | 带额外空间的径向数组容量 |
| `max_orbitals` | 127 | 轨道容量 |
| `max_orbitals_minus_one` | 126 | 轨道容量减一 |
| `max_orbitals_minus_two` | 125 | 轨道容量减二 |

在 `compile_time_capacities.f90` 修改表中第一列的基础常量；工作区容量和两个轨道容量派生量自动计算。
无需修改转导模块。修改后须全量重编译。`max_radial_points` 与运行时点数 `radial_point_count` 不同。

## 各入口的原版默认值

四个保留的 RMCDHF 入口及 RHFS、RIS、RDENSITY 各自初始化网格；
g_J 继续复用 RHFS。运行时变量仍保存在原有 `radial_grid_state`、`atomic_calculation_state` 和 `radial_node_count_settings`。

| 核模型 | RNT | H | N | HP |
| --- | --- | --- | --- | --- |
| 有限核 | `2D-6/Z` | `0.05` | `max_radial_points`（默认 590） | 0 |
| 点核 | `EXP(-65D0/16D0)/Z` | `0.0625` | `min(220,max_radial_points)` | 0 |

上表保留原版公式记号。当前代码中 `RNT/H/N/HP/Z` 分别为
`radial_grid_scale/transformed_grid_step/radial_point_count/asymptotic_grid_step/nuclear_charge`。

RMCDHF 的 `THRESH=0.05D0`、`METHOD=3`、`ODAMP=0` 保持原版。
旧 orbopt 程序及其专用实验开关已移除。

`atomic_calculation_state` 以 `working_tolerance` 保存交互提示和补丁配置中的 ACCY。初始化时序沿用原版：

- RMCDHF 在可选网格交互之前赋值 `ACCY=H**6`；交互修改 H 后不会自动重算。
  用户可在随后的 ACCY 提示中显式设置。续算应明确输入原 ACCY。
- RHFS/RIS/RDENSITY 在网格交互之后赋值 `ACCY=H**6`。
- 已移除集中初始化附带的 `UPDATE_RADIAL_ACCURACY` 和 `VALIDATE_RADIAL_GRID` 调用。
  输入设置必须满足原版的网格容量和数值约束。

输入 `.w` 经插值到本次运行网格，不会自动将 N/H/RNT 改成文件网格。
已否决的 strict/diag 独立求解器已[移除](rejected_rmcdhf_variants.md)；现行入口只覆盖保留程序。

CSFG 的专用源码、构建和测试已移除，其历史参数见[归档记录](../archive/graspg_integration.md)。

普通 `rwfnestimate`、串行/MPI `rci`、串行/MPI `rbiotransform` 和三个 `rtransition` 变体新增八个网格入口；普通 RCI 保留原有公共网格初始化，未启用注释中的点核分支。`rangular` 不初始化径向网格，但通过同一 `grasp_parameters` 使用轨道容量和打包参数。普通 `rwfnestimate` 的六个实现/接口文件已移除局部径向容量定义；修改中央容量无需逐个修改这些文件。

## 网格实验与重编译

在独立源码副本上使用 [网格脚本](../../scripts/README.md)。`atomic-suite` 布局
修改 `compile_time_capacities.f90` 和需要变更的十五个普通入口；显式 ACCY 修改也覆盖两个普通 RCI 续算入口。
仅使用 `--nnnp` 时只改中央配置的 `max_radial_points`，工作区容量自动派生，转导模块及入口保持原样；
`--n` 单独选择实际点数，其他网格参数仍修改各入口。脚本检查六项定义与完整转导接口。

容量修改示例（只预览，不写入）：

```fish
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --layout atomic-suite --grasp . --nnnp 1990 --diff
```

追加 `--apply` 才会写入；示例值用于说明接口，当前配置仍为 590。

不传数值项表示保留；`--accy 0` 在 suite 布局中写回 `H**6` 公式。

脚本显式修改 H（含点核 H）或 ACCY 时，会在 RMCDHF 的 ACCY 提示前同步阈值，
用户随后仍可显式覆盖；这是主动应用网格补丁的行为，未经补丁的源码保留原版时序。
仅修改容量不会额外更改 ACCY 时序。

恢复本次容量后，请使用新的构建目录全量编译；旧 build/bin 或已安装 bin 不会自动更新：

```fish
cmake -S . -B build-upstream-parameters -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-upstream-parameters --parallel 4
ctest --test-dir build-upstream-parameters --output-on-failure
```

使用 `build-upstream-parameters/bin/` 下的本次产物，或自行安装后核实可执行文件路径。
Linux 服务器需在服务器上重新编译。`shared_grid_capacity` 检查共享数组编译容量一致；
网格脚本回归覆盖原版布局和 suite 十五个网格入口、两个续算入口的预览、应用、重复应用、备份与恢复。
工程测试不能代替生产模型的收敛验证；候选网格仍按
[网格实验方案](../plans/grasp_grid_parameter_tuning.md) 验收。

## 参数恢复与转导验证（2026-10-13）

- 六个参数默认值与参考原版相同，通过转导模块公开；八个入口除普通 MPI 原有 PED 初始化及空白外与原版一致。
- 串行/MPI 清理后全量重编译成功；运行时共享容量检查输出 `590 600`。
- 12 项 CTest 通过，包括严格停止逻辑、RDENSITY 串行及 1/2/4 进程回归；
  完整测试在允许本地 MPI socket 通信的环境运行。
- 网格脚本 30 个离线用例通过，含真实八入口源码副本的应用与逐字节恢复，
  以及在临时副本中仅修改配置后编译验证 590/600、1990/2000 的转导和共享数组容量；
  原版容量的 `--nnnp 590 --check` 返回无待修改文件。
- 2 项严格内存 MPI 重启测试在 macOS 明确跳过：原版内存检测依赖 Linux `/proc/meminfo`。
  生产 AS6 未运行，也未提交集群作业。
