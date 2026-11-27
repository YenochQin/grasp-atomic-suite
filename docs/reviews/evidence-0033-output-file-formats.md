# Suite 生成文件的命名、结构与读写关系

2026-11-27 核对 Suite `bd4dcaa` 及其上包含状态名输入改造的工作树。范围为当前 CMake 声明的
25 个应用目标，含串行、MPI 和 OpenMP 变体。现行行为来自文件打开、写入和读取路径的静态核查。
文末的候选命名和迁移顺序供后续设计讨论。 [应用构建入口](../../CMakeLists.txt)

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

## 双正交轨道与混合系数

`rbiotransform` 和 `rbiotransform_mpi` 读取两态的描述符文件对、`.wfn` 和
`.mix` 或 `.cimix`。它们先旋转两侧径向轨道，再对各侧 CSF 混合系数施加对应的反变换。
三个跃迁程序以同一对描述符和这两套变换产物计算矩阵元。MPI 版在启动目录由 rank 0
写共享产物，角缓存和排序临时文件使用各 rank 的工作目录。
[串行入口](../../src/apps/rbiotransform/rbiotransform_main.f90) ·
[MPI 入口](../../src/apps/rbiotransform_mpi/rbiotransform_mpi_main.f90) ·
[跃迁输入](../../src/apps/rtransition/README.md)

| 现行文件 | 数据与生命周期 |
| --- | --- |
| `A.bw`、`B.bw` | 按侧保存变换后的轨道，调用公共 `RAD--WFN` 写入器。 |
| `A.bm`、`B.bm` | 以 MCDHF `.mix` 为输入，写 `CSF--MIX`。沿用输入的块目录、基组摘要、根编号和能量，负载为反变换后的系数。 |
| `A.cbm`、`B.cbm` | 以 RCI `.cimix` 为输入，使用与 `.bm` 相同的流格式。后缀与应用选项共同表达混合系数来源。 |
| `A.TB`、`B.TB` | 可复用的一体角系数。按对称块和 κ 组写 Fortran 顺序无格式记录，包含积分与系数计数、打包轨道对、系数累积终点、角系数和 bra/ket CSF 编号。尾记录 `BIOTB001` 验证描述符摘要、尺寸、κ 分组、所选 `J`、MPI 进程数和 rank。无可用缓存时重算。关闭保存时用匿名临时文件传递重算结果。 |
| `fort.<80+k>` | κ 组的排序前临时记录，含块内 bra/ket、打包轨道对和一体系数。排序后关闭并删除。 |

[变换轨道写入](../../src/apps/rbiotransform/write_transformed_radial_wavefunctions.f90) ·
[变换混合系数写入](../../src/apps/rbiotransform/read_biorthogonal_mcp.f90) ·
[角缓存记录](../../src/apps/rbiotransform/sort_biorthogonal_mcp.f90) ·
[角缓存校验](../../src/lib/core/biorthogonal_angular_cache.f90)

### `RAD--WFN` 径向流

初猜 `.initial.wfn`、优化轨道 `.wfn`、双正交轨道 `.bw` 和自然轨道 `.nw` 共用这一布局。
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

`.mix`、`.cimix`、`.bm`、`.cbm` 共用小端无记录标记的流格式。整数为有符号 int64，
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

### 双态交接的改造建议

1. **将共享产物绑定到有顺序的态对。** 在同一目录先做 `A→B` 再做 `A→C`，
   后一次会替换 `A.bw` 和同来源的 `A.bm` 或 `A.cbm`。A 的描述符摘要可以相同，
   此时它无法区分两次旋转。轨道先写，系数后写，途中失败也可能留下混合批次。
   建议一次双态输入使用一个结果目录，目录内使用 `initial.biorth.wfn`、`final.biorth.wfn`
   及对应的 `.biorth.mix` 或 `.biorth.cimix`。同一态对需并存两种来源时可使用独立结果目录。
   如需在读入时拒绝串档，可在两侧文件完整写出后发布清单，记录初末态、基组摘要、轨道和混合输入身份及对应输出。
   轨道和系数继续使用已有公共流。 [轨道写出顺序](../../src/apps/rbiotransform/rbiotransform_main.f90) ·
   [跃迁端混合读取](../../src/apps/rtransition/calculate_oscillator_strengths.f90)
2. **同态跃迁可消除磁盘副本。** 目前三个跃迁程序都会创建 `S_CP.bm` 或 `S_CP.cbm`，
   逐字节复制原流。串行及相位版正常收尾会删除，MPI 收尾没有执行删除。当前代码用独立的可分配路径
   构造副本名，23 字符状态名的截断已修复。后续可用两个独立读取句柄访问同一文件，保留两侧的块遍历位置。
   [复制程序](../../src/apps/rtransition/prepare_transition_mixing_coefficients.f90) ·
   [串行收尾](../../src/apps/rtransition/calculate_oscillator_strengths.f90) ·
   [MPI 收尾](../../src/apps/rtransition_mpi/calculate_oscillator_strengths.f90)
3. **角缓存可随其读写路径独立迁移。** 若重做 `.TB`，可按对称块和 κ 组建目录，
   使用显式字段位宽、偏移和完成标记，并保留可复用性检查。触发这一步的依据应来自跨构建复用、
   随机读取或 I/O 成本的需求。 [现行缓存键](../../src/lib/core/biorthogonal_angular_cache.f90)

## 求解器和角积分程序的其他文件

| 应用 | 现行文件 | 内容与依赖 |
| --- | --- | --- |
| `rwfnestimate` | `S.initial.wfn`、非默认设置下的 `S.rwfnestimate.summary.txt`、可选 `S.rwfnestimate.dbg` | 初猜轨道走 `RAD--WFN`，供四个 RMCDHF 入口读取。摘要和调试为文本。 [程序入口](../../src/apps/rwfnestimate/rwfnestimate_main.f90) |
| `rangular` | `mcp.30`、`mcp.31`、`mcp.(32+k)`、`S.tmp.NN`、`S.rangular.log`、可选 `S.rangular.dbg`、非默认路径的 `fort.24`、可选计时和复用统计 CSV | `mcp.30` 保存块信息和稀疏列结构，`mcp.31` 保存一电子 T 系数，`mcp.(32+k)` 保存秩 `k` 的二电子 V 系数。稀疏结构、积分组及其系数使用顺序无格式记录，排序前记录在临时通道。两个串行 RMCDHF 程序读取这套文件。 [通道头](../../src/apps/rangular/connect_mcp_files.f90) · [系数记录](../../src/apps/rangular/sort_mcp_coefficients.f90) |
| `rangular_mpi` | 各 rank 的 `mcp.mcpidx`、发布暂存文件 `mcp.mcpidx.part-<批次>`、`S.tmp<RRR>.NN`、`S.rangular.log`、可选 `S.rangular.dbg`、`fort.24` 和 `<前缀>.rank<RRR>.csv` | `mcp.mcpidx` 已有 `GMCPIDX1` 头、基组摘要、rank 布局、64 位块及组目录、稀疏索引和分组系数，尾标为 `MCPDONE1`。两个 MPI RMCDHF 程序按相同基组、批次和布局读取分片。 [MPI MCP 格式](../implemented/design-0020-rmcdhf-mcp-io.md) · [分片实现](../../src/lib/mpi/mpi_indexed_mcp.f90) |
| `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi` | `S.wfn`、OL/EOL 的 `S.mix`、`S.rmcdhf.log`、`S.rmcdhf.summary`、可选 `S.rmcdhf.dbg` | 轨道和系数分别为 `RAD--WFN` 和 `CSF--MIX`。应答、计算设置、结果摘要及调试均为文本。 [通用命名](../../src/lib/core/calculation_state_files.f90) · [SCF 链](../implemented/design-0019-rmcdhf-flow-and-state.md) |
| RMCDHF 附加输出 | `npol.summary`、`mcp.io.csv`、`mcp.bits`、`fort.96`～`fort.98` | `npol.summary` 是核极化设置与势构造摘要，由启用核极化的 `rmcdhf`、`rmcdhf_mem` 或 `rmcdhf_mem_mpi` 写出。普通 `rmcdhf_mpi` 提供各 rank 的 MCP I/O CSV 和浮点位模式文本对照。最后三项是串行稀疏矩阵调试文本。 [核极化适用范围](../implemented/design-0022-rmcdhf-nuclear-polarization.md) · [MCP 诊断](../../src/apps/rmcdhf_mpi/rmcdhf_mcp_io.f90) · [矩阵调试](../../src/apps/rmcdhf/write_hamiltonian_records.f90) |
| `rci`、`rci_mpi` | `S.cimix`、`S.csum`、`S.clog`、`rci.res` 或各 rank 的 `rci<RRR>.res` | `.cimix` 为 `CSF--MIX`，`.csum`、`.clog` 为文本摘要和交互日志。`R92RES` 重启文件使用顺序无格式记录，保存基组计数、核模型、Hamiltonian 选项、网格及轨道，以及逐块的稀疏矩阵和闭壳层 Breit 位移。 [串行重启入口](../../src/apps/rci/connect_ci_restart_file.f90) · [MPI 重启入口](../../src/apps/rci_mpi/connect_ci_restart_file.f90) · [重启读取器](../../src/lib/core/rci_restart_records.f90) |

## 跃迁、能级和属性文件

| 应用 | 现行文件 | 内容与依赖 |
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

1. **先完成双正交交接和跃迁报告命名。** 在两个 `rbiotransform` 和三个 `rtransition`
   之间统一双态路径及产物身份，处理同态副本，为含带符号矩阵元列的主报告使用独立路径。
   验证应覆盖同一初态配不同末态、交换两态、同态、MCDHF/CI 两种来源，以及部分写出后重跑。
   现有辐射跃迁和双正交回归可检查能量、系数、重叠与相位是否保持。
   [双正交回归](../../test/test_biorthogonal_states.py) · [跃迁回归](../../test/test_radiative_transitions.py)
2. **再处理自然轨道和被程序读取的文本。** 自然轨道可采用 `S.natural.wfn` 和
   `S.ci.natural.wfn` 区分来源。若需要复用自然占据数，再定义与轨道身份绑定的独立输出。
   修改 `.i` / `.ci`、`.lsj.lbl` 或跃迁报告时，需同步检查字段、单位、能级键和固定列。
   工作区内的 GraspKit 跃迁 loader 按 `.t`、`.ct`、`.t.lsj`、`.ct.lsj` 识别报告，
   `MixCoefLoader.load_suite` 有单独的文件名校验，Tools 的收集器和模板也按已知名字找产物。
   [自然轨道实测范围](evidence-0032-density-natural-orbitals.md) · [RIS4 与 FICAL](evidence-0031-isotope-shifts.md)
3. **按复用和恢复需求逐项处理缓存、重启和 MCP。** `.TB`、`±nT`、`.IOB`、`.ITB`
   和 `.lsj.T` 对应不同的负载及索引访问方式。迁移时保留每种算符、基组、截断和 rank 布局的真实约束，
   并评估为 `.lsj.T` 增加准确基组身份。`R92RES` 优先考虑固定文件名冲突，布局变更需验证
   重启时的 Hamiltonian、轨道和分片能完整恢复。串行 MCP 若换成目录流，生成器和两个串行求解器
   同步验证。MPI MCP 已有目录及身份校验，可单独调整命名。 [双正交缓存验证](evidence-0027-biorthogonal-states-and-cache.md) ·
   [RCI 重启回归](../../test/test_rci_restart.py) · [MPI MCP 回归](../../test/test_rmcdhf_mcp_paths.py)
