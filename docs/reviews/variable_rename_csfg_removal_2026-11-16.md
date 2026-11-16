# 明确候选变量重命名与 CSFG 移除

日期：2026-11-16。实现基准：`c6548d651b909849340e7b5a7457d7681192804e`，本记录描述其上的工作树改动。

范围来自[变量清单](../plans/non_csfg_variable_inventory_2026-11-16.md)。用户确认本轮只处理其中
529 项明确候选，另明确要求直接移除 CSFG 代码。其余 19,579 项保持待核实，不补造名称。

## 已完成的重命名

- **528 项重命名，1 项保留**：`def_C::pi` 已具有明确意义，沿用现名。
- 746 个普通 Fortran 文件更新了 15,620 处绑定名称或导入引用，包含实现、接口、模块成员、
  派生类型字段、过程回调及 OpenMP 变量子句。没有新增或删除原清单中的数据/过程形参实体。
- 以文件和作用域解析声明、`USE`、已有导入别名、宿主关联和类型字段；没有全仓库替换短字符串。
  既有显式 `USE` 本地别名保留，只更新其指向的远端名称。
- 字符串及普通注释保持原内容；其中的历史缩写仍可能出现。OpenMP 指令中的变量引用随实际变量更新。
  对变长后的代码行补充必要续行，没有改变运算符、数值常量、类型/kind、维度含义或初始化值。
- 普通共享状态仍由原模块持有；没有引入第二份状态、指针代理或重复库。

完整的 529 项结果、当前声明位置及原名见
[机器可读映射](../reference/non_csfg_variable_renames_2026-11-16.json)。
原清单第 6 节逐项标记已实施/保留；原行号和源码哈希仍表示采集时的历史快照。

### 五个冲突候选的调整

| 作用域 | 原名 | 原候选 | 最终名称 | 原因 |
| --- | --- | --- | --- | --- |
| `mcp` | `key` | `orbital_index_radix` | `coefficient_label_radix` | 区分本地常量与导入的同名容量参数 |
| `setcon` | `ainfcm` | `bohr_radius_cm` | `reference_bohr_radius_cm` | 区分参考常量与由其赋值的模块变量 |
| `setcon` | `ccmps` | `speed_of_light_cm_per_s` | `reference_speed_of_light_cm_per_s` | 同上 |
| `setcon` | `emeamu` | `electron_mass_u` | `reference_electron_mass_u` | 同上 |
| `setcon` | `empamu` | `proton_mass_u` | `reference_proton_mass_u` | 同上 |

调整后，作用域内候选名与声明/导入绑定的冲突预检为零。`setcon` 的这些参考数值与单位均未改变。

### 容量与工具兼容性

六个容量常量的实际定义仍只有一份，位于
[`suite_parameters_M.f90`](../../src/lib/libmod/suite_parameters_M.f90)。
[`parameter_def_M.f90`](../../src/lib/libmod/parameter_def_M.f90) 同时转导语义名称和旧名别名：
`keyorb`、`nnnp`、`nnn1`、`nnnw`、`nnnwm1`、`nnnwm2`。
这六个旧名与新名指向同一命名常量，不是重复的数值配置。

[`patch_grasp_grid.py`](../../scripts/patch_grasp_grid.py) 同时识别原版和语义名称，
保留 `--nnnp/--n/--h/...` 及 TOML 配置契约。语义源码写回语义名称；字符串不参与名称翻译。
`atomic-suite` 布局移除四个 CSFG 网格入口和两个续算入口，现覆盖十五个普通网格入口及两个普通 RCI 续算入口。
当前字段、默认值与手工修改位置见[共同参数说明](../implemented/common_parameters.md)。

其他共享模块的旧变量名及已重命名的形参关键字不再作为公开旧名提供。仓库内调用和测试已同步；
外部 Fortran 消费者如果按旧名 `USE ... ONLY` 或按旧形参名称传入关键字，需要按映射更新。
程序名、文件名、数据协议字段和交互输入契约未变。

## CSFG 移除

移除 **964 个专用源码、构建及测试文件**，其中 `src/` 下 **948 个**：

- 所有 `src/appl/*csfg*`、`src/lib/*csfg*`；
- `src/tool/graspg/` 与仅供该扩展使用的 `src/lib/libdvd90mpi/`；
- `cmake/CSFGExtension.cmake`、`test/graspg/` 和专用 F2018 测试；
- 普通测试中的 `RK_CSFG` / `DAVIDSON_MPI` 专用分支，以及 CMake 注册和验收清单中的对应目标。

`libmod`、`lib9290`、`libdvd90`、`libmcp90`、`librang90`、`mpi90` 等普通共享库保留。
当前源码树已无 CSFG 专用代码路径；构建目标及测试清单也不再包含这些程序。
`GRASP_BUILD_CSFG=ON` 明确报错。OFF、`--no-csfg` 及历史验证 JSON 的固定 `csfg=false`
维度仅为旧命令/告警预算标识保留，不能重新启用扩展。

旧[集成记录](../archive/graspg_integration.md)归档；上游来源、署名和历史哈希资料保留。
本轮没有删除本地历史构建目录、用户计算数据或其他仓库的源码。

## 验证结果

使用 GNU Fortran 15.3.0；显式设置 `FC` 和 `OMPI_FC` 为 `/usr/bin/gfortran-15`。
全新 Release 与 Check 目录均启用普通串行、MPI、jj2lsj OpenMP、rangular OpenMP。
没有安装或覆盖已有程序目录。

| 检查 | 结果 |
| --- | --- |
| 普通源码解析及声明核对 | 1,880 份文本全部解析；20,108 个实体一一对应；528 项更名、pi 保留 |
| 非目标声明核对 | 没有额外声明更名、增加或丢失 |
| 字符串/普通注释 | 746 个修改文件逐段内容一致；兼容模块新增说明注释单独登记 |
| Release 全量编译 | 通过 |
| Release CTest | 首轮 65/67；两个名称相关维护检查适配后定向重跑 2/2，最终 67 项通过，无跳过 |
| Check 全量编译及 CTest | 编译通过，67/67 通过，无跳过 |
| 网格脚本专项 | 36 个测试通过，含原版/语义名称、容量转导、真实入口修改与精确恢复 |
| CSFG 旧入口 | 临时目录配置 `GRASP_BUILD_CSFG=ON` 被预期错误拒绝 |
| 来源快照负例 | 代码、字符串、注释及前一来源 manifest 的篡改均被拒绝 |

两个维护检查分别是 OpenMP 持久缓存白名单 `raca → small_six_j_cache`，以及 MPI 接口负例
生成器中 `ierr → mpi_error_code`。它们仍检查同一缓存边界及有效/无效 ABI 调用，没有放宽断言。
十份 Fortran 测试通过 `USE` 别名适配共享变量的新公开名称，保留各测试自己的局部命名。

RMCDHF 来源校验新增
[变量重命名快照](../reference/rmcdhf_variable_rename_snapshot.json)，
逐一绑定原小写化快照的 480 份导入文件身份，再校验当前精确 SHA-256。
原上游、F2018 兼容修改和小写化快照未覆盖、未修改。其他 12 份导入文件仍沿用原检查。

### 与提交前二进制直接对照

基线是提交 `c6548d6` 对应的小写化后 Release 构建。使用既有合成测试输入，在短临时路径运行
两版二进制，以避免原有 MPI 路径长度限制；结果保存在本地审计目录。

| 用例 | 配置 | 字节一致的输出文件对 |
| --- | --- | ---: |
| rangular | 串行 | 6 |
| rangular_omp | 4 线程 | 6 |
| rangular_mpi | 2 ranks | 12 |
| rdensity | 加权多块/主量子数间隔两组 | 4 |
| rdensity_mpi | 同上，2 ranks | 4 |
| jj2lsj | even/odd/f3，4 线程 | 12 |
| 合计 | 10 组 | 44 |

MCP、密度/自然轨道、带符号转换系数和标签等输出逐字节一致。
这支持本轮变量重命名的软件行为保持，不代表完整原子 SCF/RCI 或物理收敛性验收。

### 已有验证限制

历史告警预算仍因源码覆盖缺口拒绝验收：现有编译清单含 1,911 个不同源码文件，
预算仅覆盖 1,874 个，差异为小写化阶段已经存在的 37 个文件。
本轮与小写化后的 Check 构建源码清单一致，告警数量均为 2,024，分类均为
243 条 correctness-review、1,569 条 general、212 条 obsolescent-standard。
只核对了总数与分类，不能据此声称逐条告警身份完全相同；重命名和续行也会改变诊断名称及定位。
本轮未更新历史预算，也未声称这些既有风险已解决。

`git diff --check` 仍指出普通 `rmcdhf`/`rmcdhf_mpi` 的 `in.f90` 中两处原有行尾空白；
与基准对应行核对后确认并非新增。其他修改文件未引入空白检查问题。

未运行完整 36 配置矩阵或生产规模 SCF；旁存 profiling 实现未作为独立应用构建。
本轮只处理 529 项明确候选，未核实项仍留待后续审查。

## 本地审计记录

`build-variable-rename-20261116/` 被 Git 忽略，包含本轮原文备份、转换脚本、预览、构建和数值输出：

- `source-before/`、`source-before.json`、`csfg-removal.json`：原文与删除范围；
- `proposals.json`、`rename-edits.json`、`candidate-adjustments.json`、`declaration-audit.json`：逐项映射和校验；
- `release/`、`check/` 及 `logs/`：编译、CTest/JUnit、告警与拒绝旧 CSFG 入口的日志；
- `differential-summary.json`、`differential-verified/`：44 对输出；
- `final-validation.json`：最终汇总。

本记录和映射文件随源码保存；本地构建目录、日志及数值输出不纳入版本控制。
