# grasp-atomic-suite

全部保留的 CSF 消费程序现已使用原生描述符 Parquet + 配套 header TOML。
输入约定、Rust 依赖与测试方式见[原生输入说明](docs/implemented/design-0001-native-csf-descriptors.md)。
reader 源码已内置于本仓库，随 suite 编译为 `build/lib/libcsf_descriptor.a`，
无需外部 rCSFs 源码或预编译 reader 库，计算运行不依赖 Python 转换。

中文 | [English](README.md)

`grasp-atomic-suite` 是一个面向相对论原子结构与原子性质计算的 Fortran 研究开发仓库。它保留原版 RMCDHF，并与 `rhfs`、`ris4`、`gj` 共享 GRASP 数值库，用于轨道优化、超精细结构、同位素位移和 Landé `g_J` 因子的计算、验证与性能优化。

本仓库不是上游 GRASP 的最小镜像。它保留了传统 GRASP 数值库和应用组织方式，同时加入了 CMake 构建、MPI 目标、`gj` 独立程序、回归数据和面向 `g_J`/MPI 优化的开发文档。

普通应用目录现已去掉旧 `90` 标记，Landé 因子程序命令为 `gj` 和 `gj_mpi`；对应关系和重新构建方式见[程序命名说明](docs/implemented/design-0005-program-naming.md)。

## 代码结构

```text
.
├── configure.sh              # CMake out-of-source 配置入口
├── CMakeLists.txt            # CMake 总构建文件
├── src/
│   ├── apps/
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
| 普通 GRASP 计算链 | `src/apps/rangular*`、`rwfnestimate`、`rci*`、`rbiotransform*`、`rtransition*` | `rangular{,_mpi}`、`rwfnestimate`、`rci{,_mpi}`、`rbiotransform{,_mpi}`、`rtransition{,_phase,_mpi}` | 十个目标复用 Suite 共享库与集中容量参数。[导入记录](docs/archive/2026-11-08-plan-0002-standard-grasp-integration.md) |
| `jj2lsj_omp` | `src/apps/jj2lsj` | `jj2lsj_omp` | 普通 CSF 的 OpenMP jj–LS 转换，按 LS 行并行并保留原截断顺序；见[构建、运行和验证说明](docs/implemented/design-0028-jj2lsj-openmp.md)。与 `jj2lsj_csfg` 分开。 |
| `rlevels` | `src/apps/rlevels` | 从 RMCDHF `.mix` 和 RCI `.cimix` 生成能级表，可附加 `jj2lsj` 标签。见[用法和输出格式](docs/implemented/design-0025-rlevels.md)。 |
| `gj` | `src/apps/gj` | `gj`, `gj_mpi` | 从 RHFS 计算链路中拆出的 Landé `g_J` 因子程序，读取 `isodata`、该状态的 CSF 描述符文件对、`name.mix/name.cimix`、`name.wfn`，输出 `name.gj/name.cgj`。 |
| `rhfs` | `src/apps/rhfs` | `rhfs`, `rhfs_mpi` | 相对论超精细结构程序，计算超精细常数和相关矩阵元，输出 `name.h/name.ch` 与 `name.hoffd/name.choffd`。 |
| `rdensity` | `src/apps/rdensity` | `rdensity`、`rdensity_mpi` | 径向电子密度与按 `2J+1` 平均的自然轨道，详见[迁移与验证说明](docs/implemented/design-0026-rdensity-mpi.md)。 |
| `ris4` | `src/apps/ris4` | `ris4`, `ris4_mpi` | 相对论同位素位移程序，计算正常质量位移、特殊质量位移和场位移电子因子，输出 `name.i/name.ci` 及中间角向数据。 |
| 原版 RMCDHF + 可选 NP | `src/apps/rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | 从工作区 grasp 新复制，保留原版求解流程；仅这四个读取 `npol.inp`。详见[NP 接口与来源](docs/implemented/design-0022-rmcdhf-nuclear-polarization.md)。 |

已否决的 strict/diag 独立求解器及专用测试、运行入口均已[移除](docs/implemented/adr-0002-retire-strict-and-diagnostic-rmcdhf.md)。

默认 CMake 配置启用 MPI 目标，要求配套的 MPI Fortran 工具链提供 `mpi_f08` 模块。可通过 `GRASP_ENABLE_MPI=OFF` 配置串行构建。串行轨道与性质程序需要显式开启。`jj2lsj_omp` 默认构建，`GRASP_JJ2LSJ_OPENMP=ON` 在可用时启用 OpenMP。线程数由 `OMP_NUM_THREADS` 控制，变换文件读写使用串行路径。

## 构建

服务器的 `gfortran` 已升级至 GNU Fortran 16.2.1（Red Hat 16.2.1-3），可以直接编译 Suite。
先加载配套的 MPI 和 BLAS/LAPACK 环境，再从 Suite 根目录运行以下命令，确认
`gfortran` 和 Open MPI 的 `mpifort` 都报告主版本 16。

```fish
set -gx FC (command -v gfortran)
set -gx OMPI_FC $FC
$FC -dumpfullversion
mpifort -dumpfullversion
./configure.sh --build-dir build-gf16
cmake --build build-gf16 --parallel 4
```

如果服务器使用带版本号的命令名，将 `gfortran` 换成 `gfortran-16`。
从旧工具链切换时使用新构建目录。在 GraspKit 工作区中，将构建目录放在
`../temp/tests/` 下。直接调用 CMake 时传入 `-DCMAKE_Fortran_COMPILER=$FC`。

CMake 只启用 Fortran，并用它探测 BLAS 和 LAPACK。
Cargo 使用 Rust 1.99.0 或更新版本构建内置 reader，并调用 C 编译器构建其 Zstd 依赖。
Cargo 始终使用 `--locked`，默认允许下载缺失的依赖。设置
`-DGRASP_CSF_READER_OFFLINE=ON` 可仅使用本地缓存构建。
`RCSFS_TEST_CLI` 只用于开发测试离线生成夹具，不参与 reader 构建和生产读取。
运行包含描述符夹具的 CTest 前，将支持描述符 V3 的 `rcsfs` 加入 `PATH`，
或在配置时传入 `-DRCSFS_TEST_CLI=/path/to/rcsfs`。

默认安装前缀为 `<build>/install`；可用 `--prefix PATH` 或 `-DCMAKE_INSTALL_PREFIX=PATH` 指定独立位置。
普通构建不安装，也不覆盖源码 `bin/`。

`configure.sh` 默认传入 `-DGRASP_ENABLE_MPI=ON`。如果检测到 MPI Fortran，安装后会得到 `bin/gj_mpi`、`bin/rhfs_mpi`、`bin/ris4_mpi`、新复制的 `rmcdhf_mpi`、`rmcdhf_mem_mpi`以及对应库文件。

同时构建四个原版 RMCDHF 和性质计算程序：

```fish
cmake -S . -B build-all -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
  -DCMAKE_Fortran_COMPILER=$FC -DRCSFS_TEST_CLI=/path/to/rcsfs
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
  -DCMAKE_Fortran_COMPILER=$FC -DRCSFS_TEST_CLI=/path/to/rcsfs
cmake --build build-serial --target install -j4
```

相同源码和编译器可以更新现有配置。迁移前配置了外部 reader 的缓存需要删除
`RCSFS_DESCRIPTOR_SOURCE_DIR`、`RCSFS_DESCRIPTOR_LIBRARY`，或改用新构建目录。
需要清理构建产物时运行：

```fish
./scripts/clean-build-artifacts.sh
```

### 严格 F2018 验收

沿用上面的 `FC` 和 `OMPI_FC` 设置，把支持 V3 的 `rcsfs` 加入 `PATH`，然后在
GraspKit 工作区内的 Suite 根目录配置独立的 Check 构建。服务器现有的 GNU Fortran 16.2.1
可直接使用默认告警基线。

```fish
set build_dir ../temp/tests/suite-gf16-check
./configure.sh --debug --profile f2018-check --build-dir "$build_dir" \
  -DRCSFS_TEST_CLI=(command -v rcsfs)
../graspkit-tools/.venv/bin/python test/f2018/check.py "$build_dir" --jobs 4
```

优化后的严格验收使用 `--profile f2018-audit`，省略 `--debug`，并另建构建目录。
验收会核查实际编译选项、目标和告警预算，然后运行注册的 CTest。[F2018 测试说明](test/f2018/README.md)

## 依赖与配置

构建需要 Fortran 编译器、CMake 3.31 或更新版本，以及使用 32 位 Fortran 整数的 LP64 BLAS/LAPACK。CMake 在配置时验证所选库的接口。MPI 目标还需要提供 `mpi_f08` 的 MPI Fortran 工具链。

CMake 会自动查找 BLAS/LAPACK，并要求 GNU Fortran 15 或更新版本。
严格验收默认使用 GNU Fortran 16.2 系列的告警基线。
默认 `f2018-release` 使用 O2；`f2018-check` 和 `f2018-audit` 分别提供
标准/接口/运行时检查与优化后诊断。legacy、隐式 SAVE 及参数不匹配放行选项已移除。
现行约束见[构建策略](docs/implemented/design-0002-fortran-build-profiles.md)和
[兼容性与共享接口](docs/implemented/design-0004-f2018-compatibility.md)；
当前原生输入依赖及验证边界见[描述符说明](docs/implemented/design-0001-native-csf-descriptors.md)。
历史测试数量不代表当前源码或物理精度已通过验收。

## 输入与输出

所有应用在交互读取状态名时都显示 `Name of state:`。双正交变换和跃迁程序分别标明初态与末态，FICAL 分别标明上态与下态。MPI 版本仅由 rank 0 询问。支持命令行状态名的程序可直接传入状态名。[公共状态输入实现](src/lib/core/calculation_state_files.f90)

在计算目录依次执行 `rwfnestimate`、`rangular`、`rmcdhf`，在每步的 `Name of state:` 提示输入相同状态名，例如 `Zn_odd_AS1`。第一步直接读取 `Zn_odd_AS1_header.toml`，后续步骤沿用同一前缀，计算文件在各步间直接交接。

| 步骤 | 按状态名读取 | 按状态名写出 |
| --- | --- | --- |
| `rwfnestimate` | `name_header.toml` | 初猜轨道 `name.initial.wfn` |
| `rangular` | `name_desc.parquet` 和配套头文件 | 角系数 `mcp.NN`、`name.rangular.log` |
| `rmcdhf`、`rmcdhf_mem` | 同一描述符文件对、初猜轨道和 MCP | 优化轨道 `name.wfn`、OL/EOL 混合系数 `name.mix`、`name.rmcdhf.summary`、`name.rmcdhf.log` |
| `rci` | 同一描述符文件对和 `name.wfn` | CI 混合系数 `name.cimix`、`name.csum`、`name.clog` |

MPI 版本由 rank 0 询问状态名。`rangular_mpi` 在每个 rank 的临时目录写入 `mcp.mcpidx`，随后用相同基组和 rank 布局运行 `rmcdhf_mpi` 或 `rmcdhf_mem_mpi`。已有旧格式的作业需用当前 `rangular_mpi` 重新生成。描述符和轨道位于永久 I/O 目录，角系数程序和 SCF 的日志位于启动目录。`isodata` 和适用的 `npol.inp` 按目录共用。`rci` 沿用先询问默认设置、再询问状态名的顺序。供 `gj` 等性质程序使用的状态名应控制在 23 个字符以内。[MPI MCP 格式和验收](docs/implemented/design-0020-rmcdhf-mcp-io.md)

RMCDHF 和 RCI 的串行与 MPI 版本都要求为 CSF 基组中的每个 `J/Parity` 块选择至少一个有效 ASF 根。每块输入一行根编号。空行、`0` 或其他无效选择会终止运行并指出对应的块。RMCDHF 将选中的各根纳入共同的轨道优化。选择多个根时，还需输入 EOL 权重选项。旧的仅对角 (E)AL MCP 数据需要用 Suite 对应的 `rangular` 程序重新生成角系数。

每个工作目录共用一套 MCP。应先完成本基组的 RMCDHF，再在同一目录为下一基组生成角系数。并行计算不同基组时，分别使用独立的工作目录和 MPI scratch 目录。`A` 和 `A.initial` 两个状态也应分目录保存，因为两者都会写 `A.initial.wfn`。下游程序读取 `.wfn` 轨道文件、`.mix` 的 RMCDHF 系数或 `.cimix` 的 RCI 系数。

状态名留空会使用 `rcsf_desc.parquet`、`rcsf_header.toml`、`rwfn.inp`、`rwfn.out` 和 `rmix.out`。旧输入格式以 y/n 设置回答开头。若状态名恰为单字母 Y 或 N（含小写），相应输入目录里的同名头文件优先。固定文件名的已有计算仍可通过 [save_rmcdhf 归档](scripts/README.md#保存-rmcdhf-计算结果)。

常用输入文件：

- `isodata`：核参数和同位素相关数据。
- `name_desc.parquet` 和 `name_header.toml`：CSF 描述符及其配套头文件。
- `name.wfn`：径向轨道波函数。
- `name.mix`：MCDHF/非 CI 混合系数文件。
- `name.cimix`：CI 混合系数文件。

常用输出文件：

- `gj`：`name.gj` 或 `name.cgj`。
- `rhfs`：`name.h/name.ch` 和 `name.hoffd/name.choffd`。
- `ris4`：`name.i/name.ci`，并可能生成 `name.IOB`、`name.ITB` 等角向中间文件。

`gj_mpi` 接受状态名和 `--nonci` 或 `--ci`。从计算目录运行，并将
`suite_bin` 设为构建目录下 `bin` 的绝对路径。CI 模式读取 `name.cimix`。

```fish
mpiexec -n 2 "$suite_bin/gj_mpi" name --nonci
mpiexec -n 2 "$suite_bin/gj_mpi" name --ci
```

不带参数时，程序回到交互式输入流程。

### 径向波函数二进制流

初猜轨道、优化轨道、双正交轨道和自然轨道都采用 `RAD--WFN` 二进制流。文件依次包含 16 字节文件头、每轨道 32 字节的参数表和 `float64` 径向矩阵。数值采用小端序。矩阵沿径向点连续存储，Python 形状为 `(2*K+1, L)`。第 0 行保存共享网格，后续按轨道顺序排列 `P` 和 `Q`。每条轨道仅前 `point_count` 点有效。读取器要求 `RAD--WFN` 签名。旧 `G92RWF` 顺序记录需重新生成 Suite 流文件。具体字段和 NumPy 示例见 [二进制布局](README.md#radial-wavefunction-stream)。

### ASF 混合系数二进制流

RMCDHF 的 `.mix` 与 RCI 的 `.cimix` 共用 `CSF--MIX`。文件依次包含基组 SHA-256、
对称性块目录、ASF 根编号与绝对总能量以及混合系数。每个 ASF 的有符号
`float64` 系数连续存放，列顺序与对应 V3 描述符块一致。双正交输出 `.bm` 和 `.cbm` 使用相同格式，保留来源基组的 SHA-256。
读取器要求 `CSF--MIX` 签名。旧 `G92MIX` 顺序记录需重新生成 Suite 流文件。[字段布局和 NumPy 示例](README.md#indexed-asf-mixing-stream)

## 验证数据

`data/` 目录包含名为 `test` 的示例数据集：

```text
data/isodata
data/test_desc.parquet
data/test_header.toml
data/test.mix
data/test.wfn
data/test.gj
```

复现示例 Landé 因子时，使用配套的混合文件、描述符及头文件。按实际构建目录设置 `suite_bin`。

```fish
set suite_root (pwd)
set suite_bin (realpath build-gf16/bin)
mkdir -p ../temp/tests
set gj_work (mktemp -d ../temp/tests/gj-example.XXXXXXXX)
cp data/isodata data/test_desc.parquet data/test_header.toml data/test.mix data/test.wfn "$gj_work/"
cd "$gj_work"
printf '%s\n' y test n | "$suite_bin/gj"
diff -u "$suite_root/data/test.gj" test.gj
cd "$suite_root"
```

CTest 覆盖描述符输入、计算流程、公共库、MPI 稀疏缓冲区和网格修改脚本。

```fish
ctest --test-dir build-all --output-on-failure
```

轨道回归脚本见 [test/rmcdhf](test/rmcdhf/README.md)。大型计算输入和结果留在仓库外。旧 orbopt 运行器和 SLURM 模板已移除，见[退役说明](docs/implemented/adr-0001-retire-rmcdhf-orbopt.md)。

## 径向网格与仓库迁移

公共容量集中定义在 [`cal_constants.f90`](src/lib/grasp_shared/cal_constants.f90)，
经 [`grasp_parameters.f90`](src/lib/grasp_shared/grasp_parameters.f90) 转导。库名及模块名的完整对照见[共享库命名](docs/implemented/design-0006-shared-library-naming.md)。Suite 的默认值为
`max_radial_points=590`、`radial_workspace_capacity=600`、`max_orbitals=127`、
`max_orbitals_minus_one=126`、`max_orbitals_minus_two=125` 和 `orbital_index_radix=215`。
有限核默认使用 `radial_point_count=max_radial_points`，点核分支初值为 `min(220,max_radial_points)`。
输入轨道会插值到当前计算网格。原版 RMCDHF 在网格交互之前计算 `ACCY=H**6`；
交互改变 H 时应随后明确设置 ACCY。

网格批量修改工具见 [scripts/patch_grasp_grid.py](scripts/README.md)：
复制 [scripts/config.example.toml](scripts/config.example.toml) 为同目录的 `config.toml`，
填写网格参数和运行模式后，不带参数执行即可，也可以直接传入网格参数。
本地 `config.toml` 已被 Git 忽略，省略可选参数表示不修改。
脚本根据自己的实际位置找到与 `scripts/` 同级的 `src/`。
命令行 `--nnnp` 和 TOML `nnnp` 都设置径向网格的编译期容量 `max_radial_points`。
可为不同参数保留独立源码副本，按本仓库 README 的 CMake 流程手动编译。
默认仅预览；修改后必须全量重编译，并验证物理量的网格收敛。
源码目录残留的 `.mod` 文件会覆盖新模块，CMake 会拒绝这种混合构建；
需先将旧生成文件移出源码目录。

参数列表、调用顺序和验证见[公共参数说明](docs/implemented/design-0003-shared-grid-parameters.md)。
本仓库的脚本布局按需修改集中容量配置和十五个独立网格入口。显式修改 `working_tolerance` 还覆盖
串行和 MPI 两个 RCI 续算入口。
参数调整的具体数值、原因、测试数据及验收标准见
[径向网格参数调整与验证规程](docs/plans/plan-0001-radial-grid-validation.md)。

迁移来源见 [RMCDHF 源码快照](docs/reference/rmcdhf_upstream_snapshot.json)，
保留的求解器见 [RMCDHF 流程说明](docs/implemented/design-0019-rmcdhf-flow-and-state.md)。

## 开发文档

[文档索引与实施状态审计](docs/README.md)逐份列出已实施功能、过时描述、
核查依据和待验证事项。

- 已实施说明：[公共参数](docs/implemented/design-0003-shared-grid-parameters.md)与[RHFS g_J 流程](docs/implemented/design-0027-rhfs-lande-g.md)。
- 待执行实验：[径向网格方案](docs/plans/plan-0001-radial-grid-validation.md)，物理收敛计算仍需运行。
- 已完成审查：[网格方法审查](docs/reviews/evidence-0001-radial-grid-methodology-review.md)，保留源码和诊断证据。
- 理论参考：[Landé 因子背景稿](docs/reference/design-0032-lande-g-theory.md)，算符规范尚需独立核对。
- [历史归档](docs/README.md#历史背景)，包括 gj 蓝图、轨道优化结果以及 RHFS/RIS 阶段报告。

## 开发注意事项

- 保持 Fortran 文件的本地风格，不要大规模重排旧代码格式。
- 应用源码文件和 Fortran 符号按职责使用小写 `snake_case`。迁移后的模块可保留 `*_interface` 后缀，文件按实现命名。`src/lib/grasp_shared/` 中每个文件定义一个同名模块。[仓库规范](AGENTS.md#coding-style--naming-conventions)
- 涉及数值行为的改动应附带输入数据、输出差异或回归说明。
- 修改 MPI 版本时，应同步检查串行版本和 MPI 版本的输入文件、输出文件与数值一致性。

## 许可证

本仓库采用 [MIT License](LICENSE)。

## GRASPG / CSFG 扩展

2026-11-16 已移除 CSFG 专用源码、目标和测试，普通程序使用的共享库保留。
`GRASP_BUILD_CSFG=ON` 会报错；`--no-csfg` 和 `GRASP_BUILD_CSFG=OFF` 仅作为旧构建命令的兼容写法保留。
来源见 [GRASPG 源码清单](docs/reference/graspg/source-manifest.json)。

### 来源与文献引用

CSFG 实现源自上游 [compas/graspg](https://github.com/compas/graspg)。本仓库曾对其进行迁移、
集成和适配；GRASPG 实现与 CSFG 方法归功于原作者。使用这些程序开展研究时，
请引用 GRASPG 软件论文、CSFG 方法论文以及其基础 GRASP2018 程序包：

1. Ran Si, Yanting Li, Kai Wang, Chongyang Chen, Gediminas Gaigalas, Michel Godefroid, and Per Jönsson. “Graspg – An extension to Grasp2018 based on configuration state function generators.” *Computer Physics Communications* **312**, 109604 (2025). [DOI: 10.1016/j.cpc.2025.109604](https://doi.org/10.1016/j.cpc.2025.109604).
2. Yanting Li et al. “Reducing the computational load – atomic multiconfiguration calculations based on configuration state function generators.” *Computer Physics Communications* **283**, 108562 (2023). [DOI: 10.1016/j.cpc.2022.108562](https://doi.org/10.1016/j.cpc.2022.108562).
3. C. Froese Fischer, G. Gaigalas, P. Jönsson, and J. Bieroń. “GRASP2018—A Fortran 95 version of the General Relativistic Atomic Structure Package.” *Computer Physics Communications* **237**, 184–187 (2019). [DOI: 10.1016/j.cpc.2018.10.032](https://doi.org/10.1016/j.cpc.2018.10.032).
