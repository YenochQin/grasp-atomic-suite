# Suite 生成文件的命名、结构与读写关系

2026-11-27 首次盘点以 Suite `bd4dcaa` 及其上的状态名输入改造为基准。实施复核以 Suite
`89d10d8` 及当时 Suite、GraspKit 和 Tools 工作树中的未提交改动为准。范围涵盖 CMake 声明的
25 个应用目标，包括串行、MPI 和 OpenMP 变体。标为迁移前的各节记录历史文件名，供回查旧结果。
[应用构建入口](../../CMakeLists.txt)

## 有序态对交接与跃迁报告

核查时，两个 `biorthtransform` 和三个 `rtransition` 已使用同一套路径约定。

| 项目 | 当前路径及行为 |
| --- | --- |
| 双正交结果目录 | `biorth/A/B/mcdhf/` 或 `biorth/A/B/ci/`。以程序提示中的初态、末态为序。逆向使用 `biorth/B/A/…/`，需要对该方向做双正交变换。 |
| 径向流 | `initial.biorth.wfn`、`final.biorth.wfn`，使用 `RAD--WFN`。 |
| 混合流 | `initial.biorth.mix`、`final.biorth.mix` 或 `initial.biorth.cimix`、`final.biorth.cimix`，使用 `CSF--MIX`。 |
| 完成交接 | `manifest` 记录有序态名、MCDHF/CI 来源、两侧基组摘要、四个原始文件和四个变换文件的 SHA-256。 |
| 普通和 MPI 跃迁报告 | `A.B.mcdhf.transition` 或 `A.B.ci.transition`，可选对应的 `.transition.lsj`。 |
| 相位报告 | `A.B.mcdhf.transition.phase` 或 `A.B.ci.transition.phase`，可选对应的 `.transition.phase.lsj`。 |

双正交写入器在替换产物前移除旧清单，完成四个文件后才写回清单。跃迁程序校验两侧基组和
输入、输出文件的摘要，因此需要保留原始轨道和混合流。读写两端通过结果目录内的 `.lock`
协调并发。MPI 由 rank 0 发布、核验清单并持锁，其他 rank 读取同一批产物。缺失清单或文件
不匹配时，按相同的初末态顺序和来源重跑双正交变换。

初末同态时，三个跃迁入口读取该态的原始轨道和混合流，两侧各自追踪块位置并共享按字节定位的
混合流。双正交路径和清单由受检工作树的 `src/lib/core/biorthogonal_state_files.f90` 管理。
[跃迁报告路径函数](../../src/lib/core/calculation_state_files.f90)

## 状态前缀与文件身份

`calculation_state_files` 统一处理 `Name of state:`。它先验证完整输入，再赋给应用变量。
初末态和 FICAL 的上下态均在两个名称通过校验后交付。MPI 公共封装在 rank 0 验证后广播。
旧答题文件首行 `y/n` 和空状态触发的固定文件名回退已移除，`Y`、`N` 是普通状态名。
现有缓冲容量仍决定各入口的上限。角积分、初猜、RMCDHF 和 RCI 最长 127 字符，
双正交、辐射跃迁和多数属性程序最长 23 字符，FICAL 最长 100 字符。
带位置参数的属性入口复用同一校验器。RLEVELS 还可读取显式混合文件路径或接受多个前缀。
[状态输入与路径函数](../../src/lib/core/calculation_state_files.f90) ·
[MPI 状态输入](../../src/lib/mpi/mpi_calculation_state.f90) ·
[RLEVELS 入口](../../src/apps/rlevels/rlevels.f90)

公共模块生成最常用的描述符、初猜轨道、SCF 轨道和混合系数路径。属性报告、跨态报告、
角系数缓存和重启文件仍有各自的命名规则。下文用 `S` 表示单态前缀，用 `A`、`B`
表示按程序提示顺序输入的两个态。成对列出的属性后缀以 MCDHF、CI 为序。
[状态文件表](../../src/lib/core/calculation_state_files.f90)

## 迁移前的双正交轨道与混合系数

`rbiotransform` 和 `rbiotransform_mpi` 读取两态的描述符文件对、`.wfn` 和
`.mix` 或 `.cimix`。它们先旋转两侧径向轨道，再对各侧 CSF 混合系数施加对应的反变换。
三个跃迁程序以同一对描述符和这两套变换产物计算矩阵元。MPI 版在启动目录由 rank 0
写共享产物，角缓存和排序临时文件使用各 rank 的工作目录。迁移后的入口位于
`src/apps/biorthtransform/` 和 `src/apps/biorthtransform_mpi/`。
[跃迁输入](../../src/apps/rtransition/README.md)

| 迁移前文件 | 数据与生命周期 |
| --- | --- |
| `A.bw`、`B.bw` | 按侧保存变换后的轨道，调用公共 `RAD--WFN` 写入器。 |
| `A.bm`、`B.bm` | 以 MCDHF `.mix` 为输入，写 `CSF--MIX`。沿用输入的块目录、基组摘要、根编号和能量，负载为反变换后的系数。 |
| `A.cbm`、`B.cbm` | 以 RCI `.cimix` 为输入，使用与 `.bm` 相同的流格式。后缀与应用选项共同表达混合系数来源。 |
| `A.TB`、`B.TB` | 可复用的一体角系数。按对称块和 κ 组写 Fortran 顺序无格式记录，包含积分与系数计数、打包轨道对、系数累积终点、角系数和 bra/ket CSF 编号。尾记录 `BIOTB001` 验证描述符摘要、尺寸、κ 分组、所选 `J`、MPI 进程数和 rank。无可用缓存时重算。关闭保存时用匿名临时文件传递重算结果。 |
| `fort.<80+k>` | κ 组的排序前临时记录，含块内 bra/ket、打包轨道对和一体系数。排序后关闭并删除。 |

迁移后的串行双正交目录中，`write_transformed_radial_wavefunctions.f90` 写轨道，
`read_biorthogonal_mcp.f90` 写反变换后的混合系数，`sort_biorthogonal_mcp.f90` 写排序后的角缓存。
[角缓存校验](../../src/lib/core/biorthogonal_angular_cache.f90)

### `RAD--WFN` 径向流

初猜 `.initial.wfn`、优化轨道 `.wfn`、双正交轨道和自然轨道共用这一布局。迁移前双正交轨道
使用 `.bw`，自然轨道使用 `.nw`。
数据无 Fortran 顺序记录标记，整数为小端有符号 int32，浮点数为小端 IEEE binary64。
下表的字节偏移从零开始。

| 字节偏移 | 记录 |
| --- | --- |
| 0 | 8 字节 ASCII 签名 `RAD--WFN`。 |
| 8 | 轨道数 `K`、最大有效径向点数 `L`，各为 int32。 |
| 16 | `K` 条 32 字节轨道记录。依次为 `n`、`κ`、有效点数、保留零值、轨道能量参数、原点参数。 |
| `16 + 32*K` | float64 共享径向数组。Fortran 形状为 `(L, 2*K+1)`，C/NumPy 形状为 `(2*K+1, L)`。 |

径向数组首先保存一列共享的 `r`，单位为 bohr，随后按轨道表顺序交替保存 `P` 和 `Q`。
各轨道超出其有效点数的尾部补零。轨道能量参数采用 Suite 的束缚能约定，单位为 Hartree。
文件后缀表达轨道用途，读取器按签名及记录大小识别数据布局。
[径向流读写实现](../../src/lib/core/radial_wavefunction_format.f90) ·
[格式说明及读取示例](../../README.md#radial-wavefunction-stream)

### `CSF--MIX` 混合流

`.mix`、`.cimix` 和双正交混合系数共用小端无记录标记的流格式，旧双正交后缀为 `.bm`、`.cbm`。整数为有符号 int64，
能量和有符号系数为 IEEE binary64，目录偏移和根表行号从零开始，根在原块中的编号从 1 开始。

| 区域 | 记录 |
| --- | --- |
| 112 字节文件头 | 8 字节签名 `CSF--MIX`，随后为电子数、轨道数、块数、总 CSF 数、总 ASF 根数，最后为 64 个小写 ASCII 十六进制字符表示的 V3 描述符 SHA-256。 |
| 48 字节每块 | 物理 `2J`、宇称符号、CSF 数、根数、根表首行索引、系数区的绝对字节偏移。 |
| 16 字节每根 | 原块内的根编号和绝对能量，单位为 Hartree。 |
| 系数区 | 各根的有符号系数连续存放，单根内部遵循描述符的 CSF 次序。Fortran 形状为 `(n_csfs, n_roots)`，C/NumPy 形状为 `(n_roots, n_csfs)`。 |

读取端可只访问目录和能级表。写入以 `GRMIXTMP` 标记未完成状态，完成后再写 `CSF--MIX` 签名。
描述符摘要绑定该侧的 CSF 基组，双正交变换沿用它。当前流头没有配对态、旋转批次或关联轨道的身份字段。
[混合流读写实现](../../src/lib/core/mixing_file_format.f90) ·
[格式说明及读取示例](../../README.md#indexed-asf-mixing-stream)

### 迁移前的双态交接问题

在同一目录先做 `A→B` 再做 `A→C` 时，旧写入器会替换 `A.bw` 和同来源的 `A.bm` 或 `A.cbm`。
两次旋转可以使用相同的 A 基组摘要。轨道先写、系数后写，途中失败也可能留下混合批次。
旧同态入口用 `S_CP.bm` 或 `S_CP.cbm` 逐字节复制原流。串行及相位版正常收尾会删除副本，MPI 版会保留。
[迁移后的同态读取](../../src/apps/rtransition/calculate_oscillator_strengths.f90)

## 迁移前的求解器和角积分程序文件

| 应用 | 迁移前文件 | 内容与依赖 |
| --- | --- | --- |
| `rwfnestimate` | `S.initial.wfn`、非默认设置下的 `S.rwfnestimate.summary.txt`、可选 `S.rwfnestimate.dbg` | 初猜轨道走 `RAD--WFN`，供四个 RMCDHF 入口读取。摘要和调试为文本。 [程序入口](../../src/apps/rwfnestimate/rwfnestimate_main.f90) |
| `rangular` | `mcp.30`、`mcp.31`、`mcp.(32+k)`、`S.tmp.NN`、`S.rangular.log`、可选 `S.rangular.dbg`、非默认路径的 `fort.24`、可选计时和复用统计 CSV | `mcp.30` 保存块信息和稀疏列结构，`mcp.31` 保存一电子 T 系数，`mcp.(32+k)` 保存秩 `k` 的二电子 V 系数。稀疏结构、积分组及其系数使用顺序无格式记录，排序前记录在临时通道。两个串行 RMCDHF 程序读取这套文件。 [通道头](../../src/apps/rangular/connect_mcp_files.f90) · [系数记录](../../src/apps/rangular/sort_mcp_coefficients.f90) |
| `rangular_mpi` | 各 rank 的 `mcp.mcpidx`、发布暂存文件 `mcp.mcpidx.part-<批次>`、`S.tmp<RRR>.NN`、`S.rangular.log`、可选 `S.rangular.dbg`、`fort.24` 和 `<前缀>.rank<RRR>.csv` | `mcp.mcpidx` 已有 `GMCPIDX1` 头、基组摘要、rank 布局、64 位块及组目录、稀疏索引和分组系数，尾标为 `MCPDONE1`。两个 MPI RMCDHF 程序按相同基组、批次和布局读取分片。 [MPI MCP 格式](../implemented/design-0020-rmcdhf-mcp-io.md) · [分片实现](../../src/lib/mpi/mpi_indexed_mcp.f90) |
| `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | `S.wfn`、OL/EOL 的 `S.mix`、`S.rmcdhf.log`、`S.rmcdhf.summary`、可选 `S.rmcdhf.dbg` | 轨道和系数分别为 `RAD--WFN` 和 `CSF--MIX`。应答、计算设置、结果摘要及调试均为文本。 [通用命名](../../src/lib/core/calculation_state_files.f90) · [SCF 链](../implemented/design-0019-rmcdhf-flow-and-state.md) |
| RMCDHF 附加输出 | `npol.summary`、`mcp.io.csv`、`mcp.bits`、`fort.96`～`fort.98` | `npol.summary` 是核极化设置与势构造摘要，由启用核极化的 `rmcdhf`、`rmcdhf_mem` 或 `rmcdhf_mem_mpi` 写出。普通 `rmcdhf_mpi` 提供各 rank 的 MCP I/O CSV 和浮点位模式文本对照。最后三项是串行稀疏矩阵调试文本。 [核极化适用范围](../implemented/design-0022-rmcdhf-nuclear-polarization.md) · [MCP 诊断](../../src/apps/rmcdhf_mpi/rmcdhf_mcp_io.f90) · [矩阵调试](../../src/apps/rmcdhf/write_hamiltonian_records.f90) |
| `rci`、`rci_mpi` | `S.cimix`、`S.csum`、`S.clog`、`rci.res` 或各 rank 的 `rci<RRR>.res` | `.cimix` 为 `CSF--MIX`，`.csum`、`.clog` 为文本摘要和交互日志。`R92RES` 重启文件使用顺序无格式记录，保存基组计数、核模型、Hamiltonian 选项、网格及轨道，以及逐块的稀疏矩阵和闭壳层 Breit 位移。 [串行重启入口](../../src/apps/rci/connect_ci_restart_file.f90) · [MPI 重启入口](../../src/apps/rci_mpi/connect_ci_restart_file.f90) · [重启读取器](../../src/lib/core/rci_restart_records.f90) |

## 迁移前的跃迁、能级和属性文件

| 应用 | 迁移前文件 | 内容与依赖 |
| --- | --- | --- |
| `rtransition`、`rtransition_mpi`、`rtransition_phase` | `A.B.t` / `A.B.ct`，两侧 LSJ 标签可用时另有 `A.B.t.lsj` / `A.B.ct.lsj` | 主报告按多极类型分段，记录上下能级、`J`、宇称、能差或波长、规范、跃迁率、振子强度及线强。LSJ 报告附能级标签。相位版沿用主报告名并增加带符号的矩阵元列，同目录运行普通版和相位版可能共用路径。 [报告命名](../../src/apps/rtransition/write_input_summary.f90) · [相位列](../../src/apps/rtransition_phase/README.md) |
| 三个 `rtransition` 的缓存与临时文件 | `A.B.+1T`、`A.B.-1T` 等，`fort.93`，同态的 `S_CP.bm` / `S_CP.cbm` | `±nT` 按多极算符保存顺序无格式角系数，按两侧对称块和径向积分标签分组。尾记录 `RADKT001` 绑定有序的两侧基组、算符、`J` 筛选和 rank 布局。`fort.93` 是排序前临时记录，`_CP` 是混合流副本。 [缓存实现](../../src/lib/core/transition_angular_cache.f90) · [缓存负载](../../src/apps/rtransition/sort_transition_angular_coefficients.f90) |
| `rdensity`、`rdensity_mpi` | `S.d` / `S.cd`、`S.nw`、可选 `rdensity.dbg` | 密度文本按能级分组，记录 `r`、`r²ρ(r)`、`ρ(r)`。`.nw` 是按 `2J+1` 加权平均的自然轨道，使用 `RAD--WFN`。占据数和旋转矩阵写 stdout。同名 CI 与 MCDHF 会共用 `.nw` 路径。 [密度和自然轨道写出](../../src/apps/rdensity/density_calculation.f90) |
| `gj`、`gj_mpi` | `S.gj` / `S.cgj` | 定宽文本，含核自旋、核矩，及每根的 `J`、宇称、`g_J`、`Δg_J` 和总 `g_J`。 [写入器](../../src/apps/gj/write_lande_g_factors.f90) |
| `rhfs`、`rhfs_mpi` | `S.h` / `S.ch`、`S.hoffd` / `S.choffd`、可选 `hfs92.dbg` | 对角报告含 `A`、`B` 和 Landé 因子。耦合报告含能级对的相互作用常数和按 `F` 重耦合的带符号矩阵元。 [报告命名](../../src/apps/rhfs/open_summary_file.f90) · [矩阵元写出](../../src/apps/rhfs/calculate_hyperfine_matrix_elements.f90) |
| `ris4`、`ris4_mpi` | `S.i` / `S.ci`、可选 `S.IOB`、`S.ITB` 及默认名为 `rdensity.dbg` 的调试文件 | 属性文本含 NMS、SMS 及所选场位移模型的电子因子，由 `fical` 读取。`.IOB` 和 `.ITB` 分别存一体和二体角系数的顺序无格式记录。尾标 `RISANG01` 检查基组摘要、算符、轨道标签基数和截断值。 [RIS4 流程](../../src/apps/ris4/README.md) · [缓存校验](../../src/apps/ris4/open_isotope_shift_angular_cache.f90) · [调试路径](../../src/apps/ris4/configure_debug_output.f90) |
| `fical` | `A.B.fi` / `A.B.cfi` | 用上下态 `.i` 或 `.ci` 及可用的 `.lsj.lbl` 写线位移文本，分列质量位移、场位移及合计。报告单位可选 GHz、MHz、meV。读取器核对输入节名、单位和能级身份。 [读写入口](../../src/apps/ris4/isotope_shift_transition_main.f90) |
| `jj2lsj_omp` | `S.lsj.lbl`，可选 `S.lsj.c`、`S.lsj.j`、`S.uni.lsj.lbl`、`S.uni.lsj.sum`、`S.lsj.T` 和 `jj2lsj_omp.profile.csv` | `.lbl` 存能级和主要 LSJ 组分，`.c` 存 LS-CSF 列表，`.j` 存分块的能级与 LS 展开向量，`.uni.*` 存唯一标签和分配摘要。`.lsj.T` 用 `jj2lsj` 标识的顺序无格式记录保存块信息及 `(jj 行, LS 行, 变换系数)`。标签由 RLEVELS、跃迁程序和 FICAL 消费。 [LSJ 写入与回读](../../src/apps/jj2lsj/jj2lsj_conversion.f90) · [计时输出](../../src/apps/jj2lsj/jj2lsj_profile.f90) |
| `rlevels` | stdout，Tools 外层可保存为 `S.level` 或 `S.clevel` | 按能量排序，报告块内根位置、`J`、宇称、绝对能量、相对最低能级的激发能、相邻能级间隔和可选 LSJ 标签。 [能级表写出](../../src/apps/rlevels/rlevels.f90) |

`fort.24`、`fort.734` 和 `fort.99` 是摘要、计算日志和调试单元尚未打开时的共享回退名。
源码另保留 `rangular` 旧重启分支中的 `fort.20`，当前串行和 MPI 主程序均关闭该分支。
`S_desc.parquet`、`S_header.toml`、`isodata`、`npol.inp` 和 `disks`
在这里作为计算输入统计。 [摘要与日志回退](../../src/lib/grasp_shared/calculation_file_units.f90) ·
[调试回退](../../src/lib/grasp_shared/debug_output_settings.f90) ·
[角积分入口](../../src/apps/rangular/rangular_main.f90)

## 后续迁移顺序

三项主体改动已在受检工作树接通。下列顺序沿用原核查提出的实施次序。

1. **双正交交接和跃迁报告命名。** 两个 `biorthtransform` 和三个 `rtransition`
   已统一双态路径和产物身份。同态直接读取原始文件，相位报告另有路径。回归覆盖同一初态配不同末态、
   交换两态、同态、MCDHF/CI 两种来源，以及部分写出后重跑。含点号的态对仍可能共用报告路径，见下节。
   [双正交回归](../../test/test_biorthogonal_states.py) · [跃迁回归](../../test/test_radiative_transitions.py)
2. **自然轨道和程序读取的文本。** RDENSITY、RIS4 和 JJ2LSJ 按混合系数来源写出下文所列文件。
   Suite、GraspKit 和 Tools 的相关读取、作业检查和收集路径已同步。自然轨道仍有两种路径重合情形。
   [密度回归](../../test/rdensity/regression.py) · [RIS4 与 FICAL](evidence-0031-isotope-shifts.md)
3. **角缓存、RCI 重启和 MCP。** 已按用途、状态或态对命名，保留各自的复用和恢复条件。
   JJ2LSJ 变换缓存新增 V3 基组摘要、每块保留行和完成记录校验。`R92RES`、串行 MCP 和其余角缓存
   沿用顺序记录，MPI MCP 沿用索引流 `GMCPIDX1`。受检工作树的
   `docs/implemented/design-0034-cache-restart-mcp-names.md` 记录路径及归档文件迁移规则。
   [双正交缓存验证](evidence-0027-biorthogonal-states-and-cache.md) ·
   [RCI 重启回归](../../test/test_rci_restart.py) · [MPI MCP 回归](../../test/test_rmcdhf_mcp_paths.py)

### 已复现的路径冲突

- **自然轨道可覆盖其他状态的输入或结果。** 在同一临时目录放入合法状态 `case` 和 `case.natural`
  的描述符、混合流和轨道。对 `case` 运行 RDENSITY 的 MCDHF 模式后，读回 `case.natural.wfn`，
  确认该路径原有的输入轨道被自然轨道替换。对 `case` 运行 CI 模式，再对具有不同混合系数的
  `case.ci` 运行 MCDHF 模式，也确认后一次替换了前一次写出的 `case.ci.natural.wfn`。
  串行和 MPI 入口共用自然轨道路经及写入逻辑。
  [路径函数](../../src/lib/core/calculation_state_files.f90) · [径向文件写入](../../src/apps/rdensity/density_calculation.f90)
- **含点号的不同态对可覆盖同一份跃迁报告。** 在单电子夹具中依次计算
  `left.dot → right` 的 E1 和 `left → dot.right` 的 M1。两对各有双正交清单，
  角缓存路径也含有初态名长度。两次都写到 `left.dot.right.mcdhf.transition`，
  第二份报告替换了第一份。普通、MPI、相位及对应 LSJ 报告共用将两态名以点号直接相接的路径函数。
  [报告路径函数](../../src/lib/core/calculation_state_files.f90) · [角缓存路径函数](../../src/lib/core/transition_angular_cache.f90)

### 本机验证范围

GNU Fortran 16.2.1、Open MPI 5.0.11 的 `f2018-check` 构建通过，包含 25 个应用目标。
37 项 Suite 定向用例通过，覆盖双正交及跃迁、密度和自然轨道、RIS4 与 FICAL、JJ2LSJ、
RLEVELS、角缓存、RCI 重启及串行和 MPI MCP。部分 MPI 用例触发了已有的 128 字符路径缓冲限制，
缩短夹具路径后通过。另在串行夹具中确认双正交写入器和跃迁读取器会等待冲突的锁，锁释放后完成。

GraspKit 和 Tools 共通过 556 项 Python 用例，其中 30 项使用本次构建的 Suite 校验描述符作业流程。
上述结果来自本机固定夹具和工作流回归。两处路径覆盖由单独的临时夹具复现。

## 自然轨道与文本交接的现行路径

右两列分别对应从 `.mix` 和 `.cimix` 得到的结果。

| 产物 | MCDHF | CI |
| --- | --- | --- |
| RDENSITY 自然轨道 | `S.natural.wfn` | `S.ci.natural.wfn` |
| RIS4 电子因子 | `S.mcdhf.isotope-shift` | `S.ci.isotope-shift` |
| JJ2LSJ 标签 | `S.mcdhf.lsj.labels` | `S.ci.lsj.labels` |

新写出的自然轨道沿用 `RAD--WFN`。占据数和旋转矩阵仍随 stdout 输出。
RIS4、JJ2LSJ 和跃迁文本沿用既有字段、单位和固定列。FICAL、RLEVELS 和三个跃迁入口
按使用的混合系数来源选取文本。优先读取上表中的文件名，再尝试旧版带 `.txt` 的同名文件。
两者均缺失时，RIS4 电子因子可读旧 `.i` 或 `.ci`，LSJ 标签可读旧 `.lsj.lbl`。
[Suite 路径函数](../../src/lib/core/calculation_state_files.f90)

GraspKit 的普通跃迁目录扫描支持上文列出的新文件名、带 `.txt` 的旧版文件名及旧 `.t`、`.ct`、`.t.lsj`、`.ct.lsj`。
同一态对和报告类型同时存在新旧文件时，
扫描选择新文件，也可按显式路径读取相位版。Tools 在 Suite 作业、ML 检查、能级 CSV 和
结果收集中识别新名称。上游 GRASP 作业仍使用 `.lsj.lbl`。
