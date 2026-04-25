# grasp-atomic-suite

中文 | [English](README.md)

`grasp-atomic-suite` 是一个面向相对论原子结构与原子性质计算的 Fortran 研究开发仓库。它基于 GRASP92/GRASP2018 风格代码树，当前重点维护 `rhfs90`、`ris4`、`gj90` 及其 MPI 相关组件，用于超精细结构、同位素位移和 Landé `g_J` 因子的计算、验证与性能优化。

本仓库不是上游 GRASP 的最小镜像。它保留了传统 GRASP 数值库和应用组织方式，同时加入了 CMake 构建、MPI 目标、`gj90` 独立程序、回归数据和面向 `g_J`/MPI 优化的开发文档。

## 代码结构

```text
.
├── configure.sh              # CMake out-of-source 配置入口
├── CMakeLists.txt            # CMake 总构建文件
├── src/
│   ├── appl/
│   │   ├── gj90/             # Landé g_J 独立计算程序
│   │   ├── rhfs90/           # 相对论超精细结构程序
│   │   └── ris4/             # 相对论同位素位移程序
│   └── lib/
│       ├── libmod/           # 全局参数、公共块和共享状态模块
│       ├── lib9290/          # GRASP92 风格 I/O、常数、网格和通用例程
│       ├── libdvd90/         # 对角化相关例程
│       ├── libmcp90/         # MCP 相关支持例程
│       ├── librang90/        # 角向代数、Racah/张量矩阵元例程
│       └── mpi90/            # MPI 文件、路径和并行辅助例程
├── docs/                     # 理论说明、实现追踪和优化报告
├── data/                     # gj90/RHFS 验证用示例输入与输出
├── bin/                      # 安装后的可执行文件
└── lib/                      # 安装后的静态库和 Fortran module
```

## 主要程序

| 程序 | 源码目录 | CMake 目标 | 说明 |
| --- | --- | --- | --- |
| `gj90` | `src/appl/gj90` | `gj90`, `gj90_mpi` | 从 RHFS 计算链路中拆出的 Landé `g_J` 因子程序，读取 `isodata`、`name.c`、`name.m/name.cm`、`name.w`，输出 `name.gj/name.cgj`。 |
| `rhfs90` | `src/appl/rhfs90` | `rhfs`, `rhfs_mpi` | 相对论超精细结构程序，计算超精细常数和相关矩阵元，输出 `name.h/name.ch` 与 `name.hoffd/name.choffd`。 |
| `ris4` | `src/appl/ris4` | `ris4`, `ris4_mpi` | 相对论同位素位移程序，计算正常质量位移、特殊质量位移和场位移电子因子，输出 `name.i/name.ci` 及中间角向数据。 |

默认 CMake 配置会尝试构建 MPI 版本。如果系统没有 MPI Fortran 工具链，会跳过 MPI 目标。串行应用目标默认关闭，需要在 CMake 配置时显式开启。

## 构建

推荐使用 out-of-source CMake 流程：

```sh
./configure.sh
cmake --build build --target install -j4
```

`configure.sh` 默认传入 `-DGRASP_ENABLE_MPI=ON`。如果检测到 MPI Fortran，安装后通常会得到 `bin/gj90_mpi`、`bin/rhfs_mpi`、`bin/ris4_mpi` 以及对应库文件。

调试构建：

```sh
./configure.sh --debug
cmake --build build-debug --target install -j4
```

禁用 MPI 且构建串行应用时，直接使用 CMake：

```sh
mkdir build-serial
cmake -S . -B build-serial -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-serial --target install -j4
```

重新配置前需要删除已有构建目录，或运行清理脚本：

```sh
./scripts/clean-build-artifacts.sh
```

## 依赖与配置

构建需要 Fortran 编译器、CMake、BLAS 和 LAPACK。MPI 目标还需要 MPI Fortran 编译环境，例如 `mpifort`。

CMake 会自动查找 BLAS/LAPACK，并在使用 GNU Fortran 时为较旧的 Fortran 源码增加兼容选项：

- `-fno-automatic`
- `-fallow-argument-mismatch`

## 输入与输出

这些应用沿用 GRASP 风格的状态名输入。运行时输入 `name` 后，程序会在当前目录查找相关文件。

常用输入文件：

- `isodata`：核参数和同位素相关数据。
- `name.c`：CSF 列表与耦合信息。
- `name.w`：径向轨道波函数。
- `name.m`：MCDHF/非 CI 混合系数文件。
- `name.cm`：CI 混合系数文件。

常用输出文件：

- `gj90`：`name.gj` 或 `name.cgj`。
- `rhfs90`：`name.h/name.ch` 和 `name.hoffd/name.choffd`。
- `ris4`：`name.i/name.ci`，并可能生成 `name.IOB`、`name.ITB` 等角向中间文件。

`gj90_mpi` 支持命令行方式：

```sh
mpiexec -n 2 ./bin/gj90_mpi test --nonci
mpiexec -n 2 ./bin/gj90_mpi test --ci
```

不带参数时，程序回到交互式输入流程。

## 验证数据

`data/` 目录包含名为 `test` 的示例数据集：

```text
data/isodata
data/test.c
data/test.m
data/test.w
data/test.gj
```

这些文件可用于核对 `gj90` 的 `g_J` 输出。当前 `tests/` 目录为空；如果后续加入 CTest 脚本，可在构建目录中运行：

```sh
ctest
```

## 开发文档

- `docs/lande_g.md`：Landé `g_J` 的理论定义、MCDHF/RCI 计算思路和 LS 近似对照。
- `docs/RHFS_gJ_report.md`：追踪 `rhfs90` 中 `g_J` 相关的输入读取、径向积分、角向矩阵元和 ASF 投影过程。
- `docs/gj90_blueprint.md`：`gj90` 的独立化设计、最小依赖链和核心计算拆分方案。
- `docs/RHFS_MPI_hfsgg_optimization.md`：`HFSGG_MPI` 优化范围、问题定位和效果说明。
- `docs/1.0.1-dev.1_optimization_summary.md`、`docs/1.1.1-dev.1_optimization_summary.md`：阶段性优化总结。

## 开发注意事项

- 保持 Fortran 文件的本地风格，不要大规模重排旧代码格式。
- 新的接口/helper 单元继续使用 `*_I.f90` 命名，共享模块继续使用 `*_C.f90` 或现有目录约定。
- 涉及数值行为的改动应附带输入数据、输出差异或回归说明。
- 修改 MPI 版本时，应同步检查串行版本和 MPI 版本的输入文件、输出文件与数值一致性。

## 许可证

本仓库采用 [MIT License](LICENSE)。
