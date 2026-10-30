# GRASPG 迁移、参数归属与 CMake 构建

> 状态：2026-10-30 已实施；工程构建与小型文件格式回归，不代表真实原子 SCF/RCI 的精度或收敛验证。
> 来源：GRASPG `e157d9fc61e8d2882e0b29d644b2bb011f5f745d`；suite 修改前基线 `52c0979f5b88ba0f9806b34ec6716b805304d32b`。

## 来源与迁移范围

计算源码来自原仓库的 `srcg/appl`、`srcg/lib` 和 `srcg/tool`，分别落入
`src/appl/*csfg*`、`src/lib/*csfg*` / `src/lib/libdvd90mpi` 和 `src/tool/graspg`。
保留原作者、数值实现和原 Makefile 的程序命名；不迁入 `.git`、旧 CMake 缓存、
编译产物或 `graspgtest` 的计算输出。原 GRASPG 仓库保留为未改动的来源副本。

[原始 README](../reference/graspg/UPSTREAM_README.md) 保留上游署名和引用信息；
[来源清单](../reference/graspg/source-manifest.json) 记录原文件 SHA-256、来源提交及迁移/复用位置。
清单的哈希描述导入前源码，不是经过本次改造后的内容哈希。
原 Makefile 未启用的 27 个源文件保留用于追溯，未擅自加入默认算法路径。

上游附带的 CMake 文件缺少 `csfg_memory_man`、`csfg_tv_C`、CSFG 加载等源码，
还引用不存在的工具文件，且多个目标名与脚本使用的 Makefile 名称不同。
本次以 Makefile 实际对象清单重建 CMake，而非直接使用这些旧 CMake 文件。

## 程序与库

`GRASP_BUILD_CSFG=ON` 默认构建以下串行工具，不受 `GRASP_BUILD_SERIAL_APPS` 控制；
后一个选项控制现有串行 RMCDHF 系列。

| 程序 | 功能 |
| --- | --- |
| `rcsfggenerate_csfg` | 生成 CSFG 与标签轨道表 |
| `rcsfginteract_csfg` | 按参考空间的相互作用筛选 |
| `rcsfgzerofirst_csfg` | 零阶/一阶空间整理 |
| `rwfnestimate_csfg` | 初始波函数估计 |
| `rmixaccumulate_csfg` | 混合系数累计选择 |
| `jj2lsj_csfg` | jj 到 LSJ 转换 |
| `rcsfgexpand_csfg`、`rcsfgextend_csfg`、`rcsfgsplit_csfg` | 展开、扩展、逐层拆分 |
| `rcsfgblocksplit_csfg`、`rdistHmatrix_csfg`、`rlevels_csfg` | 分块、矩阵文件重分配、能级整理 |
| `rsave_csfg` | 安装时提供的结果保存脚本 |

启用 MPI 且找到 MPI Fortran 后，另构建 `rangular_csfg_mpi`、`rmcdhf_csfg_mpi`、
`rci_csfg_mpi`、`rci_block_csfg_mpi`。关闭 MPI 时不创建这些目标，串行工具仍可完整编译。

所有目标使用 suite 的 `mod`、`9290`、`mcp90`、`rang90`、`mpiu90` 等所需基础库；
没有另建一套基础 `libmod`、`lib9290` 或 MPI 辅助库。CSFG 新增状态、CSFG I/O、
角向及矩阵元分别组成 `mod_csfg`、`9290_csfg`、`rang90_csfg`、`me_csfg`。

专用 `dvd90mpi` 保留 GRASPG 修改的 `dvdson` 及三个接口，其余源文件直接引用
suite 的 `libdvd90` 源码。原始专用实现有 MPI 归约和至少两次迭代等差异，
不能用链接基础 `dvd90` 的方式无条件替代。

`rcsfginteract_csfg` 删除了相同的基础模块与例程副本，复用 suite 库。
其新增轨道/块状态通过 `interact_*` 模块转导基础存储；不同的调试状态、高角动量项表、
筛选例程保留应用内部实现。扩展分配器复用 `csfg_memory_man`。
每个目标使用独立 Fortran module 输出目录，防止同名接口相互覆盖。
没有将原来未启用的 OpenMP 代码路径变成默认并行实现。

## 参数修改位置

### 与 suite 共用的容量

唯一配置源是 [`suite_parameters_M.f90`](../../src/lib/libmod/suite_parameters_M.f90)，
通过 `parameter_def` 兼容接口转导：`KEYORB=215`、`NNNP=590`、`NNN1=NNNP+10`、
`NNNW=127`、`NNNWM1=NNNW-1`、`NNNWM2=NNNW-2`。

已去掉波函数估计实现/接口中私有的 `590/600/120` 定义，以及相互作用工具和
矩阵重分配工具中的私有 `2000/2010/214` 配置。所有程序必须使用与库一致的容量；
因此原独立相互作用工具的 214 轨道上限现在遵循 suite 的 127 默认值。
需要更大容量时修改中央模块并重新编译全部使用方，不能仅重编译某个应用。
原来的大网格文件也不会因为迁移而自动变成 590 点；重分配工具会在读取数组前检查容量。

`NNNP` 是编译容量，`N` 是运行点数；两者不能混同。
本次保留各程序自身的网格公式和 ACCY 初始化时序。网格 patcher 的 `atomic-suite`
布局已加入四个 CSFG 网格入口和两个 RCI 续算 ACCY 入口，共 17 个网格入口。
`--nnnp` 仍只修改集中容量；显式 `--h/--n/--accy` 才改对应入口，`--accy 0`
在普通入口和 RCI 续算入口都恢复 `H**6`。详见[共同参数说明](common_parameters.md)。

### GRASPG 专有容量与默认值

统一配置源是 [`csfg_parameters_M.f90`](../../src/lib/libmod_csfg/csfg_parameters_M.f90)：

| 配置 | 默认值 | 作用 |
| --- | --- | --- |
| `CSFG_MAX_BLOCKS` | 50 | 角向/RCI/SCF 入口、展开计数与块偏移数组 |
| `CSFG_INTERACT_MAX_BLOCKS` | 20 | 相互作用筛选的附加块表；共享 `rang_Int_C` 仍有 20 项表，扩大时需同步审查 |
| `CSFG_MAX_SPLIT_LAYERS` | 50 | split/extend 工具的层表 |
| `CSFG_LABEL_LENGTH` | 1070 | 共享标签轨道字符串 |
| `CSFG_MAX_SELECTED_LEVELS` | 100 | 混合系数累计工具的选定能级表 |
| `CSFG_LEVEL_MAX_STATES` / `CSFG_LEVEL_MAX_FILES` | 20000 / 1000 | 能级工具的状态、文件容量 |
| `CSFG_LEVEL_MAX_J` | 22 | 能级打印表；与 `JFraction` 表绑定，不能只提高常量 |
| `CSFG_INCORE_MAX_ELEMENTS` | 2147483648（64 位整数） | RCI 的 `NINCOR` 矩阵元素阈值，不是用户输入的 GB 内存上限 |
| `CSFG_MATRIX_CUTOFF` / `CSFG_RESTART_CUTOFF` | `1D-10` / `1D-20` | 两类原有矩阵筛选阈值，含义不同，分别保留 |
| `CSFG_DEFAULT_*` | 保留原值 | 对角块选项、Breit n/l 限制开关与默认值、CRITE 初值 |

运行中可变的开关仍存放在 `csfg_decide_C`，交互输入照常覆盖初值；
集中配置并没有增加新的交互问题。原来数值相同但物理/算法用途不同的常数没有合并。
局部积分系数、角动量表、文件单元号、生成器的轨道枚举表与固定格式字符串不作为任意容量旋钮。
更高 n/l、任意扩容或改变阈值仍需要检查这些约束，不能由编译成功推断其物理有效性。

修复了原有 `NCFBLKTr(20)` 与 50 块入口的不一致：该数组现在使用 `CSFG_MAX_BLOCKS`。
四处 `TotCSFs_perblock` MPI 广播改用 `SIZE(...)`，避免扩容后仍发送固定 150 项。

## 构建与使用

从 suite 仓库目录执行，默认包含 CSFG；产物位于构建目录的 `bin/` 和 `lib/`：

```sh
cmake -S . -B build-csfg -DGRASP_BUILD_CSFG=ON -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-csfg --parallel 4
ctest --test-dir build-csfg --output-on-failure
cmake --install build-csfg
```

`rsave_csfg` 对输出文件名加引号并终止选项解析，支持包含空格或以 `-` 开头的名称。

安装后程序位于 suite 的 `bin/`。也可以直接使用 `build-csfg/bin/` 中的可执行文件。
串行环境使用 `-DGRASP_ENABLE_MPI=OFF`；仅构建原有程序使用 `-DGRASP_BUILD_CSFG=OFF`。
不需要原 GRASPG 仓库存在，也不需要其传统 `make-environment` 或已编译的库。
共享库、MPI 和 BLAS/LAPACK 均由 suite 的 CMake 配置解析。

## 已完成验证与边界

- GNU Fortran 16.2.1：全量串行+MPI 构建、无 MPI 全量构建，16 个 CSFG Fortran 程序均成功链接。
- MPI 配置 CTest **31/31** 通过，无 MPI 配置 **20/20** 通过；涵盖原有回归、CSFG 容量、文件格式和扩容接口检查。MPI 测试需要允许本机 socket 通信。
- `GRASP_BUILD_CSFG=OFF` 配置通过；使用临时 `DESTDIR` 验证安装，16 个 Fortran 程序及 `rsave_csfg` 均生成，未覆盖仓库现有安装目录。
- 网格 patcher 35 项回归通过，包含真实新增入口的修改/恢复；原默认配置 `--nnnp 590 --check` 无待修改文件。
- 独立编译原始 `e157d9f` 的生成器和展开工具；6 个 CSF / 2 个块的输出与迁移版逐字节一致。
  将全部候选作为参考空间时，相互作用筛选保留相同 CSF 文件。
  另外独立编译原始相互作用工具：以 2 个 CSF 为参考，从 6 个候选中保留 5 个，迁移版结果逐字节一致。
- 在临时源码副本中编译 `NNNP=790`、`NNN1=800`、`NNNW=149`、`CSFG_MAX_BLOCKS=64`，
  检查共享数组和波函数显式接口尺寸一致；实际仓库默认值未改变。
- 新增回归位于 [`test/graspg`](../../test/graspg/)，测试只使用临时目录。

本次没有执行真实原子的 RMCDHF/RCI、完整 GRASPG 示例、性能基准或集群作业。
编译中仍有上游旧 MPI/BLAS 调用风格产生的参数类型/秩警告，沿用 suite 的兼容编译选项。
下一步的科学验收应使用同一输入、同一容量和相同库版本，对照原始 GRASPG 的
轨道、能量、矩阵维数和重启结果；本次小型文件格式检查不能替代该验收。
