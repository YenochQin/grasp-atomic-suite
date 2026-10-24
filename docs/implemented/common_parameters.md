# 参数转导层与原版默认值

> 文档状态：现行实现；2026-10-13 保留参数转导层，使用原版默认值与网格初始化。
> 参考源码：工作区 GRASP `488e875d1a7ec631baf4ad33682798bc795b7435`。
> 2026-10-24 新复制四个原版 RMCDHF 应用，现有十三个网格入口；来源见 [NP 实现说明](rmcdhf_nuclear_polarization.md)。
> 容量由集中配置定义，经 parameter_def 转导；历史实验数值不代表当前默认值。

## 参数归属

[`suite_parameters_M.f90`](../../src/lib/libmod/suite_parameters_M.f90) 是容量常量的唯一来源；
[`parameter_def_M.f90`](../../src/lib/libmod/parameter_def_M.f90) 用 `USE suite_parameters, ONLY: ...`
转导六个原有名称。共享库和各串行/MPI 程序继续 `USE parameter_def`，
修改配置不需要逐个修改调用方。配置模块默认值与参考原版相同。
网格的 `radial_grid_defaults` 集中初始化模块不使用，各入口保留原版初始化。

| 常量 | 默认值 | 含义 |
| --- | --- | --- |
| `KEYORB` | 215 | 轨道整数打包参数 |
| `NNNP` | 590 | 径向数组编译容量；部分旧径向路径也将它用于尾部边界 |
| `NNN1` | 600 | 带额外空间的径向数组容量 |
| `NNNW` | 127 | 轨道容量 |
| `NNNWM1` | 126 | `NNNW-1` |
| `NNNWM2` | 125 | `NNNW-2` |

在 `suite_parameters_M.f90` 修改基础常量。`NNN1=NNNP+10`、`NNNWM1=NNNW-1`、
`NNNWM2=NNNW-2` 自动派生；无需修改转导模块。修改后须全量重编译。`NNNP` 和运行时点数 `N` 不是同一变量。

## 各入口的原版默认值

十个 RMCDHF 入口（四个新复制的原版、四个既有 orbopt、独立严格内存 MPI 和 P0 诊断版本）及 RHFS、RIS、RDENSITY 各自初始化网格；
g_J 继续复用 RHFS。运行时变量仍保存在原有 `grid_C`、`def_C` 和 `coun_C`。

| 核模型 | RNT | H | N | HP |
| --- | --- | --- | --- | --- |
| 有限核 | `2D-6/Z` | `0.05` | `NNNP`（默认 590） | 0 |
| 点核 | `EXP(-65D0/16D0)/Z` | `0.0625` | `MIN(220,NNNP)` | 0 |

RMCDHF 的 `THRESH=0.05D0`、`METHOD=3`、`ODAMP=0` 保持原版。
普通 orbopt MPI 原有的 `PED` 初始化和算法开关保留。

ACCY 初始化时序也恢复原版：

- RMCDHF 在可选网格交互之前赋值 `ACCY=H**6`；交互修改 H 后不会自动重算。
  用户可在随后的 ACCY 提示中显式设置。续算应明确输入原 ACCY。
- RHFS/RIS/RDENSITY 在网格交互之后赋值 `ACCY=H**6`。
- 已移除集中初始化附带的 `UPDATE_RADIAL_ACCURACY` 和 `VALIDATE_RADIAL_GRID` 调用。
  输入设置必须满足原版的网格容量和数值约束。

输入 `.w` 经插值到本次运行网格，不会自动将 N/H/RNT 改成文件网格。
独立 [`rmcdhf_mem_strict_mpi`](rmcdhf_mem_strict_mpi.md) 使用原版参数初始化，
仅其停止策略及诊断输出另有改动。`GRASP_COUNT_ACCY` 等原有 orbopt 实验覆盖
仍按对应程序实现处理；严格内存版本不读取这些实验开关。

## 网格实验与重编译

在独立源码副本上使用 [网格脚本](../../scripts/README.md)。`atomic-suite` 布局
修改 `suite_parameters_M.f90` 和需要变更的十三个入口。
仅使用 `--nnnp` 时只改中央配置的 NNNP，NNN1 自动派生，转导模块及入口保持原样；
`--n` 单独选择实际点数，其他网格参数仍修改各入口。脚本检查六项定义与完整转导接口。

容量修改示例（只预览，不写入）：

```sh
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --layout atomic-suite --grasp . --nnnp 1990 --diff
```

追加 `--apply` 才会写入；示例值用于说明接口，当前配置仍为 590。

不传数值项表示保留；`--accy 0` 在 suite 布局中写回 `H**6` 公式。

脚本显式修改 H（含点核 H）或 ACCY 时，会在 RMCDHF 的 ACCY 提示前同步阈值，
用户随后仍可显式覆盖；这是主动应用网格补丁的行为，未经补丁的源码保留原版时序。
仅修改容量不会额外更改 ACCY 时序。

恢复本次容量后，请使用新的构建目录全量编译；旧 build/bin 或已安装 bin 不会自动更新：

```sh
cmake -S . -B build-upstream-parameters -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-upstream-parameters --parallel 4
ctest --test-dir build-upstream-parameters --output-on-failure
```

使用 `build-upstream-parameters/bin/` 下的本次产物，或自行安装后核实可执行文件路径。
Linux 服务器需在服务器上重新编译。`shared_grid_capacity` 检查共享数组编译容量一致；
网格脚本回归覆盖原版布局和 suite 十三入口布局的预览、应用、重复应用、备份与恢复。
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
