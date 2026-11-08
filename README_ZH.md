# grasp-atomic-suite

中文 | [English](README.md)

`grasp-atomic-suite` 是一个面向相对论原子结构与原子性质计算的 Fortran 研究开发仓库。它保留原版 RMCDHF 和独立 strict/diag 求解器，并与 `rhfs90`、`ris4`、`gj90` 共享 GRASP 数值库，用于轨道优化、超精细结构、同位素位移和 Landé `g_J` 因子的计算、验证与性能优化。

本仓库不是上游 GRASP 的最小镜像。它保留了传统 GRASP 数值库和应用组织方式，同时加入了 CMake 构建、MPI 目标、`gj90` 独立程序、回归数据和面向 `g_J`/MPI 优化的开发文档。

## 代码结构

```text
.
├── configure.sh              # CMake out-of-source 配置入口
├── CMakeLists.txt            # CMake 总构建文件
├── src/
│   ├── appl/
│   │   ├── gj90/             # Landé g_J 独立计算程序
│   │   ├── rdensity/         # 串行/MPI 径向电子密度与自然轨道
│   │   ├── rhfs90/           # 相对论超精细结构程序
│   │   ├── rmcdhf{,_mem,_mpi,_mem_mpi}/ # 从 grasp 新复制，支持可选 NP
│   │   ├── rmcdhf90_mem_{strict,diag}_mpi/ # 独立严格停止/诊断版本
│   │   └── ris4/             # 相对论同位素位移程序
│   └── lib/
│       ├── libmod/           # 全局参数、公共块和共享状态模块
│       ├── lib9290/          # GRASP92 风格 I/O、常数、网格和通用例程
│       ├── libdvd90/         # 对角化相关例程
│       ├── libmcp90/         # MCP 相关支持例程
│       ├── librang90/        # 角向代数、Racah/张量矩阵元例程
│       └── mpi90/            # MPI 文件、路径和并行辅助例程
├── docs/                     # 状态索引与分类文档
│   ├── implemented/          # 已实施功能的现行说明与代码追踪
│   ├── plans/                # 待执行的实验和验证方案
│   ├── reviews/              # 已完成审查与证据记录
│   ├── reference/            # 理论背景参考稿
│   └── archive/              # 被替代的蓝图和历史报告
├── data/                     # gj90/RHFS 验证用示例输入与输出
├── test/                     # 轨道优化和公共库测试
├── scripts/                  # 网格修改及构建清理工具
├── bin/                      # 安装后的可执行文件
└── lib/                      # 安装后的静态库和 Fortran module
```

## 主要程序

| 程序 | 源码目录 | CMake 目标 | 说明 |
| --- | --- | --- | --- |
| 普通 GRASP 计算链 | `src/appl/rangular90*`、`rwfnestimate90`、`rci90*`、`rbiotransform90*`、`rtransition90*` | `rangular{,_mpi}`、`rwfnestimate`、`rci{,_mpi}`、`rbiotransform{,_mpi}`、`rtransition{,_phase,_mpi}` | 十个目标复用 suite 共享库与集中容量参数；见[来源、构建与验证范围](docs/implemented/standard_grasp_integration.md)。 |
| `jj2lsj_omp` | `src/appl/jj2lsj90` | `jj2lsj_omp` | 普通 CSF 的 OpenMP jj–LS 转换，按 LS 行并行并保留原截断顺序；见[构建、运行和验证说明](docs/implemented/jj2lsj_openmp.md)。与 `jj2lsj_csfg` 分开。 |
| `gj90` | `src/appl/gj90` | `gj90`, `gj90_mpi` | 从 RHFS 计算链路中拆出的 Landé `g_J` 因子程序，读取 `isodata`、`name.c`、`name.m/name.cm`、`name.w`，输出 `name.gj/name.cgj`。 |
| `rhfs90` | `src/appl/rhfs90` | `rhfs`, `rhfs_mpi` | 相对论超精细结构程序，计算超精细常数和相关矩阵元，输出 `name.h/name.ch` 与 `name.hoffd/name.choffd`。 |
| `rdensity` | `src/appl/rdensity` | `rdensity`、`rdensity_mpi` | 径向电子密度与按 `2J+1` 平均的自然轨道，详见[迁移与验证说明](docs/implemented/rdensity_mpi.md)。 |
| `ris4` | `src/appl/ris4` | `ris4`, `ris4_mpi` | 相对论同位素位移程序，计算正常质量位移、特殊质量位移和场位移电子因子，输出 `name.i/name.ci` 及中间角向数据。 |
| 原版 RMCDHF + 可选 NP | `src/appl/rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | 从工作区 grasp 新复制，保留原版求解流程；仅这四个读取 `npol.inp`。详见[NP 接口与来源](docs/implemented/rmcdhf_nuclear_polarization.md)。 |
| 严格内存 SCF | `src/appl/rmcdhf90_mem_strict_mpi` | `rmcdhf_mem_strict_mpi` | 从原版独立复制的内存 MPI 程序；EOL 要求轨道与能量连续两轮同时达标，并保存逐轮波函数和收敛 CSV。详见[重启与验证说明](docs/implemented/rmcdhf_mem_strict_mpi.md)。 |
| P0 分阶段诊断 | `src/appl/rmcdhf90_mem_diag_mpi` | `rmcdhf_mem_diag_mpi` | 从原版独立复制，记录匹配、阻尼和正交化前后的系数与径向值；保留原停止判据。详见[诊断说明](docs/implemented/rmcdhf_mem_diag_mpi.md)。 |

默认 CMake 配置会尝试构建 MPI 版本。如果系统没有 MPI Fortran 工具链，会跳过 MPI 目标。串行轨道与性质程序需要显式开启。`jj2lsj_omp` 默认构建，`GRASP_JJ2LSJ_OPENMP=ON` 在可用时启用 OpenMP，线程数由 `OMP_NUM_THREADS` 控制；变换文件读写模式保留串行路径。

## 构建

推荐使用 out-of-source CMake 流程（GNU Fortran 15）：

```sh
export FC="$(command -v gfortran-15)"
export OMPI_FC="$FC"
./configure.sh
cmake --build build --parallel 4
```

默认安装前缀为 `<build>/install`；可用 `--prefix PATH` 或 `-DCMAKE_INSTALL_PREFIX=PATH` 指定独立位置。
普通构建不安装，也不覆盖源码 `bin/`。严格验收入口见 [F2018 测试说明](test/f2018/README.md)。

`configure.sh` 默认传入 `-DGRASP_ENABLE_MPI=ON`。如果检测到 MPI Fortran，安装后会得到 `bin/gj90_mpi`、`bin/rhfs_mpi`、`bin/ris4_mpi`、新复制的 `rmcdhf_mpi`、`rmcdhf_mem_mpi`以及对应库文件。

同时构建新复制的四个原版 RMCDHF、独立 strict/diag 版本和性质计算程序：

```sh
cmake -S . -B build-all -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-all --parallel 4
ctest --test-dir build-all --output-on-failure
cmake --install build-all
```

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

构建需要 Fortran 编译器、CMake 3.31 或更新版本、BLAS 和 LAPACK。MPI 目标还需要 MPI Fortran 编译环境，例如 `mpifort`。

CMake 会自动查找 BLAS/LAPACK，默认 `f2018-audit` profile 使用以下标准与兼容选项：

- `-std=f2018`
- `-fno-automatic`

当前配置已通过 gfortran 15.3.0 的全量 F2018 编译/链接，见
[兼容性修复与验证](docs/implemented/f2018_compatibility.md)。
`-fallow-argument-mismatch` 已移除，MPI/数值库的接口及数组实参已修正，见
[接口验证记录](docs/reviews/argument_interfaces_2026-11-08.md)。audit/legacy 仍保留静态存储；
36/47 是 2026-11-08 Linux 的历史结果。2026-11-09 macOS 严格 Check/Release 全量构建通过，
各 64/68 测试通过、4 项 Linux 专用测试跳过，完整验收仍失败。同 O0/O2 的 rmcdhf 一轮
全输出一致；额外 MPI Davidson 早收敛探针尚未解决，见[本轮验证](docs/reviews/f2018_strict_interfaces_npol_2026-11-09.md)。

## 输入与输出

性质计算程序沿用 GRASP 风格的状态名输入。运行时输入 `name` 后，程序会在当前目录查找相关文件。轨道优化程序沿用 RMCDHF 的交互输入，使用 `rcsf.inp`、`rwfn.inp` 和 MCP 数据。`rwfnestimate`、`rangular` 和 `rangular_mpi` 已在本仓库构建；`rnucleus`、普通 `rcsfgenerate` 仍来自外部 GRASP，或使用已核对的输入文件。

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

这些文件可用于核对 `gj90` 的 `g_J` 输出。CTest 已接入公共库积分、MPI 稀疏缓冲区、轨道事务/轮次逻辑和网格修改脚本测试：

```sh
ctest --test-dir build-all --output-on-failure
```

轨道回归脚本见 [test/rmcdhf](test/rmcdhf/README.md)。大型计算输入和结果留在仓库外。旧 orbopt 运行器和 SLURM 模板已移除，见[退役说明](docs/implemented/orbopt_retirement.md)。

## 径向网格与仓库迁移

公共容量集中定义在 [`suite_parameters_M.f90`](src/lib/libmod/suite_parameters_M.f90)，
经 [`parameter_def_M.f90`](src/lib/libmod/parameter_def_M.f90) 转导，默认值与原版一致：
`NNNP=590`、`NNN1=600`、`NNNW=127`、`NNNWM1=126`、`NNNWM2=125`、`KEYORB=215`。
各入口保留自身的原版网格初始化：有限核默认 `N=NNNP=590`；提供点核分支的程序使用
`N=MIN(220,NNNP)`。两个 CSFG RCI 程序沿用公共网格初始化，没有独立点核分支。
输入轨道会插值到当前计算网格。原版 RMCDHF 在网格交互之前计算 `ACCY=H**6`；
交互改变 H 时应随后明确设置 ACCY。

网格批量修改工具见 [scripts/patch_grasp_grid.py](scripts/README.md)：
复制 [scripts/config.example.toml](scripts/config.example.toml) 为同目录的 `config.toml`，
修改配置中的路径、数值和运行模式后，不带参数执行即可；原有命令行接口仍保留。
本地 `config.toml` 已被 Git 忽略，省略可选参数表示不修改。
配置留在 `scripts/` 时，本仓库设置 `GRASP_SOURCE=".."`、
`SOURCE_LAYOUT="atomic-suite"`，完整原版设置 `"grasp2018"` 和对应源码路径。
可为不同参数保留独立源码副本，按原版 README 的 CMake 流程手动编译。
默认仅预览；修改后必须全量重编译，并验证物理量的网格收敛。
源码目录残留的 `.mod` 文件会覆盖新模块，CMake 会拒绝这种混合构建；
需先将旧生成文件移出源码目录。

参数列表、调用顺序和验证见[公共参数说明](docs/implemented/common_parameters.md)。
本仓库的脚本布局按需修改集中容量配置和十七个独立网格入口；显式 ACCY 修改还覆盖
两个 CSFG RCI 续算入口，保留转导接口。
参数调整的具体数值、原因、测试数据及验收标准见
[径向网格参数调整与验证规程](docs/plans/grasp_grid_parameter_tuning.md)。

迁移来源、公共库取舍和当时验证见[历史迁移记录](docs/archive/rmcdhf_migration.md)，
其中早期网格默认值和测试数量已被公共参数说明替代。
后续轨道优化开发在本仓库进行；原 `rmcdhf_test` 工作目录保留作为参考。

## 开发文档

[文档索引与实施状态审计](docs/README.md)逐份列出已实施功能、过时描述、
核查依据和待验证事项。

- 已实施说明：[公共参数](docs/implemented/common_parameters.md)与[RHFS g_J 流程](docs/implemented/RHFS_gJ_report.md)。
- 待执行实验：[径向网格方案](docs/plans/grasp_grid_parameter_tuning.md)，物理收敛计算仍需运行。
- 已完成审查：[网格方法审查](docs/reviews/grasp_grid_methodology_review.md)，保留源码和诊断证据。
- 理论参考：[Landé 因子背景稿](docs/reference/lande_g.md)，算符规范尚需独立核对。
- 历史归档：[逐份目录](docs/README.md#逐份文档结论)，包括 gj90 蓝图、迁移快照、RHFS/RIS 阶段报告；分别标明保留代码、已删除尝试和旧测量。

## 开发注意事项

- 保持 Fortran 文件的本地风格，不要大规模重排旧代码格式。
- 新的接口/helper 单元继续使用 `*_I.f90` 命名，共享模块继续使用 `*_C.f90` 或现有目录约定。
- 涉及数值行为的改动应附带输入数据、输出差异或回归说明。
- 修改 MPI 版本时，应同步检查串行版本和 MPI 版本的输入文件、输出文件与数值一致性。

## 许可证

本仓库采用 [MIT License](LICENSE)。

## GRASPG / CSFG 扩展

默认启用 `GRASP_BUILD_CSFG=ON`，构建 12 个串行 Fortran 工具、4 个 MPI 计算程序（需要 MPI），
并安装 `rsave_csfg`。这些目标使用本仓库基础库；公共容量统一在 `suite_parameters_M.f90`，
CSFG 专有容量与初值在 `csfg_parameters_M.f90`。程序清单、迁移分析、构建和验证范围见
[GRASPG 集成说明](docs/implemented/graspg_integration.md)。

`-DGRASP_BUILD_CSFG=OFF` 可关闭整个扩展；`-DGRASP_ENABLE_MPI=OFF` 保留串行 CSFG 工具。

### 来源与文献引用

CSFG 实现源自上游 [compas/graspg](https://github.com/compas/graspg)。本仓库对其进行迁移、
集成和适配；GRASPG 实现与 CSFG 方法归功于原作者。使用这些程序开展研究时，
请引用 GRASPG 软件论文、CSFG 方法论文以及其基础 GRASP2018 程序包：

1. Ran Si, Yanting Li, Kai Wang, Chongyang Chen, Gediminas Gaigalas, Michel Godefroid, and Per Jönsson. “Graspg – An extension to Grasp2018 based on configuration state function generators.” *Computer Physics Communications* **312**, 109604 (2025). [DOI: 10.1016/j.cpc.2025.109604](https://doi.org/10.1016/j.cpc.2025.109604).
2. Yanting Li et al. “Reducing the computational load – atomic multiconfiguration calculations based on configuration state function generators.” *Computer Physics Communications* **283**, 108562 (2023). [DOI: 10.1016/j.cpc.2022.108562](https://doi.org/10.1016/j.cpc.2022.108562).
3. C. Froese Fischer, G. Gaigalas, P. Jönsson, and J. Bieroń. “GRASP2018—A Fortran 95 version of the General Relativistic Atomic Structure Package.” *Computer Physics Communications* **237**, 184–187 (2019). [DOI: 10.1016/j.cpc.2018.10.032](https://doi.org/10.1016/j.cpc.2018.10.032).
