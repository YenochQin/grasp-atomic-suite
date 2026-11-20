# grasp-atomic-suite

全部保留的 CSF 消费程序现已使用原生描述符 Parquet + 配套 header TOML。
输入约定、Rust 依赖与测试方式见[原生输入说明](docs/implemented/native_csf_descriptors.md)。
reader 源码已内置于本仓库，随 suite 编译为 `build/lib/libcsf_descriptor.a`，
无需外部 rCSFs 源码或预编译 reader 库，计算运行不依赖 Python 转换。

中文 | [English](README.md)

`grasp-atomic-suite` 是一个面向相对论原子结构与原子性质计算的 Fortran 研究开发仓库。它保留原版 RMCDHF，并与 `rhfs`、`ris4`、`gj` 共享 GRASP 数值库，用于轨道优化、超精细结构、同位素位移和 Landé `g_J` 因子的计算、验证与性能优化。

本仓库不是上游 GRASP 的最小镜像。它保留了传统 GRASP 数值库和应用组织方式，同时加入了 CMake 构建、MPI 目标、`gj` 独立程序、回归数据和面向 `g_J`/MPI 优化的开发文档。

普通应用目录现已去掉旧 `90` 标记，Landé 因子程序命令为 `gj` 和 `gj_mpi`；对应关系和重新构建方式见[程序命名说明](docs/implemented/program_naming.md)。

## 代码结构

```text
.
├── configure.sh              # CMake out-of-source 配置入口
├── CMakeLists.txt            # CMake 总构建文件
├── src/
│   ├── appl/
│   │   ├── gj/               # Landé g_J 独立计算程序
│   │   ├── rdensity/         # 串行/MPI 径向电子密度与自然轨道
│   │   ├── rhfs/             # 相对论超精细结构程序
│   │   ├── rmcdhf{,_mem,_mpi,_mem_mpi}/ # 从 grasp 新复制，支持可选 NP
│   │   └── ris4/             # 相对论同位素位移程序
│   └── lib/
│       ├── grasp_shared/     # 共享状态、容量和基础接口
│       ├── core/             # CSF 输入、径向数值计算与通用例程（libcore.a）
│       ├── davidson_eigensolver/ # Davidson 本征求解库（libdavidson_eigensolver.a）
│       ├── breit_angular_support/ # Breit 角系数与外部链接入口
│       ├── spin_angular/        # 角向代数、Racah/张量矩阵元例程
│       └── mpi/              # MPI 文件、路径和并行辅助例程
├── docs/                     # 状态索引与分类文档
│   ├── implemented/          # 已实施功能的现行说明与代码追踪
│   ├── plans/                # 待执行的实验和验证方案
│   ├── reviews/              # 已完成审查与证据记录
│   ├── reference/            # 理论背景参考稿
│   └── archive/              # 被替代的蓝图和历史报告
├── data/                     # gj/RHFS 验证用示例输入与输出
├── test/                     # 轨道优化和公共库测试
├── scripts/                  # 网格修改及构建清理工具
├── bin/                      # 安装后的可执行文件
└── lib/                      # 安装后的静态库和 Fortran module
```

## 主要程序

| 程序 | 源码目录 | CMake 目标 | 说明 |
| --- | --- | --- | --- |
| 普通 GRASP 计算链 | `src/appl/rangular*`、`rwfnestimate`、`rci*`、`rbiotransform*`、`rtransition*` | `rangular{,_mpi}`、`rwfnestimate`、`rci{,_mpi}`、`rbiotransform{,_mpi}`、`rtransition{,_phase,_mpi}` | 十个目标复用 suite 共享库与集中容量参数；见[来源、构建与验证范围](docs/implemented/standard_grasp_integration.md)。 |
| `jj2lsj_omp` | `src/appl/jj2lsj` | `jj2lsj_omp` | 普通 CSF 的 OpenMP jj–LS 转换，按 LS 行并行并保留原截断顺序；见[构建、运行和验证说明](docs/implemented/jj2lsj_openmp.md)。与 `jj2lsj_csfg` 分开。 |
| `rlevels` | `src/appl/rlevels` | `rlevels` | 从 RMCDHF `.m` 和 RCI `.cm` 生成能级表，可附加 `jj2lsj` 标签。见[用法和输出格式](docs/implemented/rlevels.md)。 |
| `gj` | `src/appl/gj` | `gj`, `gj_mpi` | 从 RHFS 计算链路中拆出的 Landé `g_J` 因子程序，读取 `isodata`、该状态的 CSF 描述符文件对、`name.mix/name.cimix`、`name.wfn`，输出 `name.gj/name.cgj`。 |
| `rhfs` | `src/appl/rhfs` | `rhfs`, `rhfs_mpi` | 相对论超精细结构程序，计算超精细常数和相关矩阵元，输出 `name.h/name.ch` 与 `name.hoffd/name.choffd`。 |
| `rdensity` | `src/appl/rdensity` | `rdensity`、`rdensity_mpi` | 径向电子密度与按 `2J+1` 平均的自然轨道，详见[迁移与验证说明](docs/implemented/rdensity_mpi.md)。 |
| `ris4` | `src/appl/ris4` | `ris4`, `ris4_mpi` | 相对论同位素位移程序，计算正常质量位移、特殊质量位移和场位移电子因子，输出 `name.i/name.ci` 及中间角向数据。 |
| 原版 RMCDHF + 可选 NP | `src/appl/rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | 从工作区 grasp 新复制，保留原版求解流程；仅这四个读取 `npol.inp`。详见[NP 接口与来源](docs/implemented/rmcdhf_nuclear_polarization.md)。 |

已否决的 strict/diag 独立求解器及专用测试、运行入口均已[移除](docs/implemented/rejected_rmcdhf_variants.md)。

默认 CMake 配置会尝试构建 MPI 版本。如果系统没有 MPI Fortran 工具链，会跳过 MPI 目标。串行轨道与性质程序需要显式开启。`jj2lsj_omp` 默认构建，`GRASP_JJ2LSJ_OPENMP=ON` 在可用时启用 OpenMP，线程数由 `OMP_NUM_THREADS` 控制；变换文件读写模式保留串行路径。

## 构建

推荐使用 out-of-source CMake 流程（GNU Fortran 15）：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC $FC
./configure.sh
cmake --build build --parallel 4
```

内置 reader 需要 Rust 1.99.0 或更新的发行版、Cargo 和 C 编译器。
Cargo 始终使用 `--locked`，默认允许下载缺失的依赖。设置
`-DGRASP_CSF_READER_OFFLINE=ON` 可仅使用本地缓存构建。
`RCSFS_TEST_CLI` 只用于开发测试离线生成夹具，不参与 reader 构建和生产读取。

默认安装前缀为 `<build>/install`；可用 `--prefix PATH` 或 `-DCMAKE_INSTALL_PREFIX=PATH` 指定独立位置。
普通构建不安装，也不覆盖源码 `bin/`。严格验收入口见 [F2018 测试说明](test/f2018/README.md)。

`configure.sh` 默认传入 `-DGRASP_ENABLE_MPI=ON`。如果检测到 MPI Fortran，安装后会得到 `bin/gj_mpi`、`bin/rhfs_mpi`、`bin/ris4_mpi`、新复制的 `rmcdhf_mpi`、`rmcdhf_mem_mpi`以及对应库文件。

同时构建四个原版 RMCDHF 和性质计算程序：

```fish
cmake -S . -B build-all -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
  -DRCSFS_TEST_CLI=/path/to/rcsfs
cmake --build build-all --parallel 4
ctest --test-dir build-all --output-on-failure
cmake --install build-all
```

调试构建：

```fish
./configure.sh --debug
cmake --build build-debug --target install -j4
```

禁用 MPI 且构建串行应用时，直接使用 CMake：

```fish
mkdir build-serial
cmake -S . -B build-serial -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON \
  -DRCSFS_TEST_CLI=/path/to/rcsfs
cmake --build build-serial --target install -j4
```

相同源码和编译器可以更新现有配置。迁移前配置了外部 reader 的缓存需要删除
`RCSFS_DESCRIPTOR_SOURCE_DIR`、`RCSFS_DESCRIPTOR_LIBRARY`，或改用新构建目录。
需要清理构建产物时运行：

```fish
./scripts/clean-build-artifacts.sh
```

## 依赖与配置

构建需要 Fortran 编译器、CMake 3.31 或更新版本、BLAS 和 LAPACK。MPI 目标还需要 MPI Fortran 编译环境，例如 `mpifort`。

CMake 会自动查找 BLAS/LAPACK，并要求 GNU Fortran 15 或更新版本。
正式验收基线使用 GNU Fortran 15.3.0。
默认 `f2018-release` 使用 O2；`f2018-check` 和 `f2018-audit` 分别提供
标准/接口/运行时检查与优化后诊断。legacy、隐式 SAVE 及参数不匹配放行选项已移除。
现行约束见[构建策略](docs/implemented/build_profiles.md)和
[兼容性与共享接口](docs/implemented/f2018_compatibility.md)；
当前原生输入依赖及验证边界见[描述符说明](docs/implemented/native_csf_descriptors.md)。
历史测试数量不代表当前源码或物理精度已通过验收。

## 输入与输出

所有应用在交互读取状态名时都显示 `Name of state:`。双正交变换和跃迁程序分别标明初态与末态，FICAL 分别标明上态与下态。MPI 版本仅由 rank 0 询问。支持命令行状态名的程序可直接传入状态名。[公共状态输入实现](src/lib/core/calculation_state_files.f90)

在计算目录依次执行 `rwfnestimate`、`rangular`、`rmcdhf`，在每步的 `Name of state:` 提示输入相同状态名，例如 `Zn_odd_AS1`。第一步直接读取 `Zn_odd_AS1_header.toml`，后续步骤沿用同一前缀，计算文件在各步间直接交接。

| 步骤 | 按状态名读取 | 按状态名写出 |
| --- | --- | --- |
| `rwfnestimate` | `name_header.toml` | 初猜轨道 `name.initial.wfn` |
| `rangular` | `name_descriptors.parquet` 和配套头文件 | 角系数 `mcp.NN`、`name.rangular.log` |
| `rmcdhf`、`rmcdhf_mem` | 同一描述符文件对、初猜轨道和 MCP | 优化轨道 `name.wfn`、OL/EOL 混合系数 `name.mix`、`name.rmcdhf.summary`、`name.rmcdhf.log` |
| `rci` | 同一描述符文件对和 `name.wfn` | CI 混合系数 `name.cimix`、`name.csum`、`name.clog` |

MPI 版本由 rank 0 询问状态名。各 rank 在各自的临时目录交接 `mcp.NN`，描述符和轨道位于永久 I/O 目录，角系数与 SCF 日志位于启动目录。`isodata` 和适用的 `npol.inp` 仍按目录共用。`rci` 沿用先询问默认设置、再询问状态名的顺序。后续若要使用 `gj` 等性质程序，状态名应控制在 23 个字符以内。

每个工作目录共用一套 `mcp.NN`。应先完成本基组的 RMCDHF，再在同一目录为下一基组生成角系数。并行计算不同基组时，分别使用独立的工作目录和 MPI scratch 目录。`A` 和 `A.initial` 两个状态也应分目录保存，因为两者都会写 `A.initial.wfn`。下游程序优先读取新后缀，对应新文件缺失时可读取归档中的 `.w`、`.m`、`.cm`。

状态名留空会使用传统的 `rcsf_*`、`rwfn.inp`、`rwfn.out` 和 `rmix.out`。以 y/n 回答设置开头的旧输入脚本仍可使用。若状态名恰为单字母 Y 或 N（含小写），相应输入目录里的同名头文件优先。同名重算时，(E)AL 不更新该状态已有的混合系数文件。固定文件名的已有计算仍可通过 [save_rmcdhf 归档](scripts/README.md#保存-rmcdhf-计算结果)。

常用输入文件：

- `isodata`：核参数和同位素相关数据。
- `name_descriptors.parquet` 和 `name_header.toml`：CSF 描述符及其配套头文件。
- `name.wfn`：径向轨道波函数。
- `name.mix`：MCDHF/非 CI 混合系数文件。
- `name.cimix`：CI 混合系数文件。

常用输出文件：

- `gj`：`name.gj` 或 `name.cgj`。
- `rhfs`：`name.h/name.ch` 和 `name.hoffd/name.choffd`。
- `ris4`：`name.i/name.ci`，并可能生成 `name.IOB`、`name.ITB` 等角向中间文件。

`gj_mpi` 支持命令行方式：

```fish
mpiexec -n 2 ./bin/gj_mpi test --nonci
mpiexec -n 2 ./bin/gj_mpi test --ci
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

这些文件可用于核对 `gj` 的 `g_J` 输出。CTest 已接入公共库积分、MPI 稀疏缓冲区、轨道事务/轮次逻辑和网格修改脚本测试：

```fish
ctest --test-dir build-all --output-on-failure
```

轨道回归脚本见 [test/rmcdhf](test/rmcdhf/README.md)。大型计算输入和结果留在仓库外。旧 orbopt 运行器和 SLURM 模板已移除，见[退役说明](docs/implemented/orbopt_retirement.md)。

## 径向网格与仓库迁移

公共容量集中定义在 [`compile_time_capacities.f90`](src/lib/grasp_shared/compile_time_capacities.f90)，
经 [`grasp_parameters.f90`](src/lib/grasp_shared/grasp_parameters.f90) 转导。库名及模块名的完整对照见[共享库命名](docs/implemented/shared_library_naming.md)。默认值与原版一致：
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
- 历史归档：[逐份目录](docs/README.md#逐份文档结论)，包括 gj 蓝图、迁移快照、RHFS/RIS 阶段报告；分别标明保留代码、已删除尝试和旧测量。

## 开发注意事项

- 保持 Fortran 文件的本地风格，不要大规模重排旧代码格式。
- 新的接口/helper 单元继续使用 `*_I.f90` 命名，共享模块继续使用 `*_C.f90` 或现有目录约定。
- 涉及数值行为的改动应附带输入数据、输出差异或回归说明。
- 修改 MPI 版本时，应同步检查串行版本和 MPI 版本的输入文件、输出文件与数值一致性。

## 许可证

本仓库采用 [MIT License](LICENSE)。

## GRASPG / CSFG 扩展

2026-11-16 已移除 CSFG 专用源码、目标和测试，普通程序使用的共享库保留。
`GRASP_BUILD_CSFG=ON` 会报错；`--no-csfg` 和 `GRASP_BUILD_CSFG=OFF` 仅作为旧构建命令的兼容写法保留。
来源与旧集成过程见[归档记录](docs/archive/graspg_integration.md)，其中的启用方式已不适用于当前源码。

### 来源与文献引用

CSFG 实现源自上游 [compas/graspg](https://github.com/compas/graspg)。本仓库曾对其进行迁移、
集成和适配；GRASPG 实现与 CSFG 方法归功于原作者。使用这些程序开展研究时，
请引用 GRASPG 软件论文、CSFG 方法论文以及其基础 GRASP2018 程序包：

1. Ran Si, Yanting Li, Kai Wang, Chongyang Chen, Gediminas Gaigalas, Michel Godefroid, and Per Jönsson. “Graspg – An extension to Grasp2018 based on configuration state function generators.” *Computer Physics Communications* **312**, 109604 (2025). [DOI: 10.1016/j.cpc.2025.109604](https://doi.org/10.1016/j.cpc.2025.109604).
2. Yanting Li et al. “Reducing the computational load – atomic multiconfiguration calculations based on configuration state function generators.” *Computer Physics Communications* **283**, 108562 (2023). [DOI: 10.1016/j.cpc.2022.108562](https://doi.org/10.1016/j.cpc.2022.108562).
3. C. Froese Fischer, G. Gaigalas, P. Jönsson, and J. Bieroń. “GRASP2018—A Fortran 95 version of the General Relativistic Atomic Structure Package.” *Computer Physics Communications* **237**, 184–187 (2019). [DOI: 10.1016/j.cpc.2018.10.032](https://doi.org/10.1016/j.cpc.2018.10.032).
