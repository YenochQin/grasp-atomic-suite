# 共享库过程模块化验收

日期：2026-11-16。修改前基准：`867583d4b338d92b44ed9e3be44f8213411983e1`。
范围：`src/lib` 中全部 188 份 `_I.f90` 所声明的 200 个具体过程。
本报告记录接口与软件行为回归，不作为新原子体系的物理收敛结论。

## 迁移清单与接口

[逐过程清单](library_module_inventory_2026-11-16.csv)列出公开过程、原模块和接口文件、
现在的实现文件、库目标、模块导入位置及消费者构建目标。导入位置是保守的源码清单：
无 `only` 的导入可能未使用全部过程，同名应用自有模块须结合构建目标判断。

| 库目标 | 迁移过程数 |
| --- | ---: |
| `9290` | 74 |
| `dvd90` | 11 |
| `mcp90` | 2 |
| `rang90` | 89 |
| `mpiu90` | 24 |
| 合计 | 200 |

- 每个实现置于模块 `contains` 中，默认 `private`，显式公开原过程。
  保留原 `_i` 模块 API；实现文件保留原路径。全部旧 `_I.f90` 已删除。
- `dvdson.f90` 原含九个外部过程。现在按已有接口模块分别编译，新增
  `adds`、`dvdrvr`、`initdvd`、`mgs_nrm`、`multbc`、`newvec`、`ovflow`、`tstsel`
  实现文件。算法、局部存储、初始化及调用次序不变。
- MPI 原接口成为 `mpiu_i` 实现模块。`startmpi`、`mpix_bytes`、`stopmpi`
  放入无套件状态依赖的 `mpi_runtime_i`，由 `mpiu_i` 重导出。
  `mpi_c`、`cpath`、`sys_getwd` 直接使用这个底层模块，消除启停/路径/状态之间的环。
- GNU Fortran 15 对与 intrinsic 同名的模块过程 `count`、`pack` 的过程指针解析失败，
  用小型独立复现确认。两个实现分别命名为 `count_impl`、`pack_impl`，由原模块重导出为
  原公开名称；不增加外部转发过程或重复形参声明。新增 `f2018_library_modules`
  回归覆盖过程指针关联和关键字调用。
- 保留 `gdvd_i` 的 `dvd_operator` 抽象回调接口；`rkco_gg` 继续使用
  `cord_callback`。没有需要迁移的命名泛型接口。抽象接口描述回调契约，不是实现的副本。

对照原接口和实现检查了形参顺序、名称、类型/kind、形状、属性及函数返回类型。
字符声明的旧写法与 `character(len=...)` 等价，不批量改写。
将 `diaga3` 的实现形参 `rec` 改为原公开名称 `recc`，将 `nine` 的 `inn` 改为 `in`，
以保留关键字调用。给 `items`、`dvdson`、`gdvd`、`nine`、`sixj4`、`sixj`、
`spicmvmpi` 补回原接口的 `intent`。

**两处原接口错误不能机械保留：**

| 过程 | 原接口 | 原实现及迁移后 |
| --- | --- | --- |
| `itjj` | `kg1: intent(in)` | `intent(out)`；过程执行 `kg1=2*kg` |
| `cpath` | 三个路径参数均为 `intent(in)` | 均为 `intent(out)`；过程填充路径 |

没有改变这两处的赋值行为。生成旧接口与新模块过程的过程指针关联编译检查，
198 个过程匹配；其余两个仅报告上表已确认的 `intent` 差异。
另外逐过程比较可执行源码，扣除已记录的形参/实现名替换以及声明和 `use` 调整后，
200 个过程的计算语句一致。

## 调用与构建

- 五个库的 CMake 清单移除旧接口源，Davidson 清单增加拆出的实现，MPI 清单增加
  `mpi_runtime.f90`。插值测试只编译 `interp.f90`，不再重复编译已删除的接口源。
- `rdensity/density_parallel_mpi.f90` 的三份手写 MPI 过程接口改用 `mpiu_i`。
- `rangular`、`rangular_mpi`、`rci`、`rci_mpi`、RIS4 的 `smsnew`/`smsnew_mpi`
  将 `external cord` 改为 `use cord_i, only: cord`，继续将该过程作为 `rkco_gg` 的实参。
  `gdvd` 的矩阵乘回调及 MPI `spicmvmpi` 的实际过程参数由全量编译和 Davidson 测试覆盖。
- `sys_getwd` 的 `stopmpi(what=..., myid=...)` 关键字调用保持有效。
- 保留测试替身作用：`gdsummpi`、`itjpo`、`slater`、`dpbdt`、`quad`
  改成同 API 的测试专用模块过程；需要时拆出独立编译单元以避免测试依赖环。
- jj2lsj 的线程内核仍从共享实现源编译为独立对象，使用目标独立的模块目录；
  OpenMP 开/关两种配置均覆盖。不复制公共库或全局状态。

检查了库外的应用、测试，以及工作区可见构建文件中的套件库链接引用。
当前可见消费者全部改为模块调用，没有必须保留的原外部符号兼容入口，也没有本批排除过程。
这次迁移改变 Fortran 链接符号；不可见的第三方二进制不在验证范围，依赖这些库的外部项目
需要使用新模块重新编译和链接，不能直接替换旧静态库而保持二进制 ABI。

## 验证与证据

使用 GNU Fortran 15，显式设置 Open MPI 的 Fortran 后端。当前机器没有仓库约定的
`/usr/bin/gfortran-15`，使用同主版本的本机安装；没有改动仓库工具链规则。
所有配置均从新的源外构建目录开始。Release 使用 `-O2`；Check 使用仓库的 F2018
严格接口诊断、`-fcheck=all`、回溯及浮点异常检查。

| 配置 | 全量构建 | CTest |
| --- | --- | --- |
| 修改前 Release，串行 + MPI + jj2lsj OpenMP | 通过 | 76/76 |
| 修改后 Release，串行 + MPI + jj2lsj OpenMP | 通过 | 77/77 |
| 修改后 Check，串行 + MPI + jj2lsj OpenMP | 通过 | 77/77 |
| 修改后 Release，MPI 关闭、jj2lsj OpenMP 关闭 | 通过 | 51/51 |

新增过程 API 测试前，修改后 Release/Check 均为 76/76，关闭 MPI/OpenMP 为 50/50。
初次受限环境运行的 MPI socket 错误、macOS 临时路径别名导致的维护测试失败，
在允许本地 socket 并使用解析后的临时目录重跑后消失；基线也使用相同环境。

直接使用既有测试输入分别运行修改前后 Release 二进制，执行零容差逐字节比较：

| 输出 | 一致的文件/输出对 |
| --- | ---: |
| 串行 rangular MCP | 6 |
| 两 rank rangular MPI MCP | 12 |
| rdensity 加权多块/间隔轨道，串行及两 rank MPI 的密度和自然轨道 | 8 |
| jj2lsj even/odd/f3，四线程的标签、CSF 和带符号转换系数 | 9 |
| QUAD、局部稀疏 MPI 核、Davidson、NMTEJJ、jj2lsj 系数测试标准输出 | 5 |
| 合计 | 40 |

最后一项 jj2lsj 系数测试成功时标准输出为空；其数值由测试内部断言检查。
不声称所有可执行文件逐字节相同，也不将有限算例回归等同于所有物理体系的证明。

本次本地证据使用 `suite-module-baseline-*`、`suite-module-candidate-*`、
`suite-module-check-*`、`suite-module-serial-*` 构建及日志前缀；
`smod-final/summary.json` 记录最终输出大小与 SHA-256，目录内保留输入、输出和运行日志。
这些构建、机器配置和计算输出未加入仓库。
