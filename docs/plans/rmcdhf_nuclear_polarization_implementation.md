# rmcdhf 核极化有效势：代码实施与代码测试计划

> 状态：2026-10-24 已完成本文范围内的代码实施与代码测试；未进行物理计算验证。
> 现行接口及测试记录见 [已实施说明](../implemented/rmcdhf_nuclear_polarization.md)。
> 实现版本为首次引入现行说明的 Git 提交，查询方式见该说明开头。
> 编写日期：2026-10-24；同日按用户更正，将接入基线改为从本地 grasp 新复制的原版 rmcdhf 系列。
> 既有 rmcdhf_orbopt 系列保留原状，已撤回此前对其进行的 NP 接入。
> 源码基准：`7f9bf2d25ad1a2c6074c61c798cc686899677da4`，编写前工作树干净。
> 范围：代码接口、数据流、构建、单元测试和受控组件集成测试。
> 本文不安排真实原子/离子的 SCF 运行、能级/同位素位移对照、网格/活性空间收敛、物理基准或生产作业。

## 1. 实施目标与支持范围

在已有核 Coulomb 势初始化后，将局域标量 E1/E2 核极化势纳入共享 `NPOT_C.ZZ`，使轨道方程和一电子哈密顿量使用同一个外势。NP 默认关闭；缺少配置文件时保持既有输入与计算路径。

先从 `grasp/src/appl/rmcdhf90*` 新复制四个原版应用到下表目录，再接入 NP；不复制公共库。原版快照见[来源清单](../reference/rmcdhf_upstream_snapshot.json)。

| 源码目录 | 构建目标 | 初始化入口 |
| --- | --- | --- |
| `src/appl/rmcdhf/` | `rmcdhf` | `GETSCD` |
| `src/appl/rmcdhf_mem/` | `rmcdhf_mem` | `GETSCD` |
| `src/appl/rmcdhf_mpi/` | `rmcdhf_mpi` | `GETSCDmpi` |
| `src/appl/rmcdhf_mem_mpi/` | `rmcdhf_mem_mpi` | `GETSCDmpi` |

既有 `rmcdhf_orbopt*`、独立 strict/diagnostic 求解器、性质程序和外部 RCI 不在此次接入范围。文档明确列出上述支持目标，不宣称其他程序会读取新配置。公共 `NUCPOT` 不自动读取 NP 文件，避免改变未接入程序的语义。

以下章节保留实施前的设计措辞；所列接口现已实现，实际约束及验证结果以[现行说明](../implemented/rmcdhf_nuclear_polarization.md)为准。

## 2. 当前源码中的接入依据

| 现有位置 | 已核实行为 | 实施约束 |
| --- | --- | --- |
| [getscd.f90](../../src/appl/rmcdhf/getscd.f90) / [getscdmpi.f90](../../src/appl/rmcdhf_mpi/getscdmpi.f90) | `RADGRD` 后调用 `NUCPOT`，再调用 `SETRWFA` / `SETRWFMPI` | 新初始化放在 `NUCPOT` 后、轨道读取前；四份入口分别接入 |
| [npot_C.f90](../../src/lib/libmod/npot_C.f90) / [nucpot.f90](../../src/lib/lib9290/nucpot.f90) | `ZZ` 容量为 `NNNP`，其物理表示为 $-rV(r)$ | 不把 `VNP` 直接加到 `ZZ`；保留 `NNUC` 的原有含义 |
| [ypot.f90](../../src/appl/rmcdhf/ypot.f90) | 串行以 `YP(:N)=ZZ(:N)` 初始化 | 读取合成后的完整外势 |
| [dacon.f90](../../src/appl/rmcdhf/dacon.f90) | 非对角一电子项中读取 `ZZ(I)/C` | 不能只改 `YPOT` 而遗漏 `DACON` |
| [rinti.f90](../../src/lib/lib9290/rinti.f90) / [setham.f90](../../src/appl/rmcdhf/setham.f90) | `MODE=0` 包含外势；`SETHAM` 对角与非对角项均使用该模式 | 通过原积分链纳入 NP，不再额外加第二份 NP 矩阵；保留 `MODE=1` 动能语义 |
| [MPI ypot.f90](../../src/appl/rmcdhf_mpi/ypot.f90) / [cofpotmpi.f90](../../src/appl/rmcdhf_mpi/cofpotmpi.f90) | `YP=ZZ/NPROCS`，随后对 `YP/XP/XQ` 求和归约 | 每个 rank 保存完整 `ZZ`，只沿用既有局部贡献分配 |
| [rscfmpivu.f90](../../src/appl/rmcdhf_mpi/rscfmpivu.f90) / [cpath.f90](../../src/lib/mpi90/cpath.f90) | 串行输入来自 `permdir`，工作目录可切换到各 rank 临时目录 | NP 配置路径必须显式从 `permdir` 构造，不能依赖初始化时的当前目录 |
| [mpiu.f90](../../src/lib/mpi90/mpiu.f90) | `stopmpi` 内含 `MPI_Barrier`、`MPI_Finalize` 和普通 `stop` | 不从 root 单独调用它处理配置错误；需要确保一致退出且返回失败状态 |
| [test/CMakeLists.txt](../../test/CMakeLists.txt) | 同时注册纯代码测试与依赖外部原子数据的回归 | 新测试必须独立标记与选择，不能以运行整个 CTest 集合作为本文验收命令 |

## 3. 输入与数据契约

### 3.1 独立可选配置

采用新文件 `npol.inp`，使用 Fortran namelist `&npol ... /`。不修改 `isodata`、已有交互 stdin 顺序或旧输入脚本。

| 字段 | 约定 |
| --- | --- |
| `schema_version` | 第一版使用 `1`，文件存在时必须显式给出；不接受未知版本 |
| `enabled` | 默认 `.false.`；只有显式 `.true.` 才启用 |
| `use_e1`, `use_e2` | 独立通道开关，默认 `.false.`；启用 NP 时至少选择一个通道 |
| `units` | 启用时必须明确为 `atomic`；第一版只接收已换算的内部单位 |
| `k1`, `b1` | E1 强度与截断半径，启用 E1 时必填 |
| `k2`, `b2` | E2 强度与截断半径，启用 E2 时必填 |
| `isotope_z`, `isotope_a` | 启用时必填；核电荷须与已读入的运行核模型一致，质量数不从平均原子量猜测 |
| `model_id` | 启用时固定为 `flambaum2021_scalar_e1_e2`；第一版拒绝其他模型标识 |
| `parameter_source` | 启用时必填；记录参数来源与 E2 约定，不触发在线查找或参数拟合 |

最小关闭示例：

```fortran
&npol
  schema_version = 1
  enabled = .false.
/
```

解析要求：

- 文件不存在视为关闭；文件存在但无法读取、语法错误、未知字段或未知版本必须报错，不能伪装成关闭。
- 开启的通道要求强度有限且非负、截断半径有限且严格为正；显式零强度合法，缺失值不能用零默认值代替。用独立的“已赋值”状态或未设置哨兵区分二者。
- 未启用通道不要求其数值，且对总势严格贡献零。整体关闭时不使用通道数值；仍检查文件语法及版本。
- 核标识和字符串长度、`units`、数值类型均有确定的错误消息。拒绝静默路径或字符串截断。
- 第一版不在程序内接收多套 fm/MeV/B(EL) 输入并自动推导模型。输入生产者负责按既定 E1/E2 约定换算，程序负责检查与记录单位。

### 3.2 内部算符表示

代码内使用

$$
V_1(r)=-\frac{K_1}{r^4+b_1^4},\qquad
V_2(r)=-\frac{K_2}{r^6+b_2^6},\qquad
V_{\mathrm{NP}}=V_1+V_2,
$$

其中 $r,b_1,b_2$ 以 $a_0$ 表示，$K_1$、$K_2$ 分别以 $E_h a_0^4$、$E_h a_0^6$ 表示。`k1/k2` 已包含文献势中 $e^2\alpha/2$ 的单位换算，不在内核里再乘一遍该系数。E2 的修正极化率与截断参数必须由输入来源保持配套。

合成规则固定为

$$
ZZ_{\mathrm{total}}(i)=ZZ_{\mathrm{Coul}}(i)-R(i)V_{\mathrm{NP}}(i).
$$

- 保存纯 Coulomb 快照 `ZZ_COUL`，另存 `VNP_E1`、`VNP_E2`、`VNP_TOTAL`；它们使用公共 `parameter_def` 容量，活动范围为 `1:N`。
- 关闭或所有已选通道均显式为零强度时，使用直接复制路径保留原 `ZZ`，不通过额外浮点运算改动默认路径。
- `R(1)=0` 时直接求有限截断势，并用 `-R*VNP` 合成；不计算 `ZZ/R`。
- `N` 越界、负半径、非有限网格/结果必须明确失败；不得访问 `NNNP` 之外或写入 `N+1:` 的数组尾部。
- 核电荷 `Z`、Fermi 参数、`NNUC`、原网格和 `ACCY` 初始化时序保持原义。

## 4. 代码拆分与修改清单

### 4.1 新增模块

| 拟新增文件 | 职责 |
| --- | --- |
| `src/lib/libmod/npol_C.f90` | 配置/状态类型、显式默认关闭状态、共享势数组或其受控访问接口；不执行文件 I/O 或 MPI |
| `src/lib/lib9290/npol_input.f90` | 读取 namelist，检查必填项、模型标识与单位，返回 `status/message`；不在库层直接 `STOP` |
| `src/lib/lib9290/npol_kernel.f90` | 纯 E1/E2 求值与按已保存基线合成总势；独立于轨道、CSF 和 MPI |
| `src/lib/lib9290/npol_init.f90` | 串行初始化编排、Coulomb 快照管理、启用状态与输出元数据 |
| `src/lib/mpi90/npol_init_mpi.f90` | root 文件 I/O、配置/势数组广播、错误同步；复用串行解析和内核 |

优先使用模块过程提供显式接口；若采用现有外部例程风格，则同时提供对应 `_I.f90`，不引入隐式过程调用。MPI 依赖只能进入 `mpiu90` 适配层，不能反向加入 `libmod` 或 `lib9290`。

`NUCPOT`、`NPOT_C` 的既有对外接口保持兼容。NP 配置只由接入的应用入口显式加载；新模块被链接到其他应用并不自动启用 NP。

### 4.2 状态生命周期

提供可分别测试的操作，建议名为：

1. `NPOL_RESET`：清理上一轮状态，恢复关闭与零数组。
2. `NPOL_READ_CONFIG`：仅解析并返回配置及错误。
3. `NPOL_VALIDATE_CONFIG`：检查完整性、核标识和单位。
4. `NPOL_BUILD_GRID`：只依据配置与显式传入网格产生 E1/E2 数组。
5. `NPOL_COMPOSE`：从显式传入或已保存的纯 Coulomb 数组生成总 `ZZ`。
6. `NPOL_INIT` / `NPOL_INIT_MPI`：应用可调用的组合入口。

捕获基线和合成必须分开：只有刚执行完 `NUCPOT` 才捕获 `ZZ_COUL`；再次 `COMPOSE` 仍读取这份快照，不能把已加 NP 的 `ZZ` 当成新基线。新核模型或网格初始化时先 reset，再重新生成和捕获基线。

内核对关闭通道先短路，不对未提供的截断参数进行除法。对半径幂采用可检查的求值策略，避免中间幂溢出产生未报告的 NaN/Inf；不改变既定势的数学定义。候选配置和数组先在临时状态中完成构建，全部成功后才发布；失败时保留本轮纯 Coulomb `ZZ`，NP 状态保持关闭，不留下部分总势。所有已验证配置与数组在 SCF 中保持只读，轨道变化不触发重新叠加。

### 4.3 四个程序入口

四个应用分别修改：

- `rscfvu.f90` / `rscfmpivu.f90`：构造明确的 NP 文件路径。
- `getscd.f90` / `getscdmpi.f90`：在 `NUCPOT` 后调用新初始化。
- `getscd_I.f90` / `getscdmpi_I.f90`：若增添 `NPFILE` 参数，同步更新接口与全部对应调用点。

建议为四个 `GETSCD*` 接口新增显式 `NPFILE` 参数。串行传入本次输入目录的 `npol.inp`；MPI 由主程序用 `permdir` 构造。所有 rank 的字符变量先初始化，非 root 不读取该路径。不要通过修改共享 `NUCPOT` 让全部程序隐式读取文件。

串行顺序：

```text
SETISO / 参数输入
SETQIC -> RADGRD -> NUCPOT
NPOL_INIT(NPFILE, ...)       # 拟新增
SETRWFA
原 MATRIX / NEWCO / SCF 路径
```

MPI 顺序相同，但轨道读取例程为 `SETRWFMPI`，新增入口为 `NPOL_INIT_MPI`。新路径必须在首次矩阵构建及 `COFPOT` 调用前完成。

### 4.4 保留既有消费者与变分逻辑

- `YPOT`、`DACON`、`RINTI(...,0)` 继续读取完整 `ZZ`。通常无需在这些消费者中再插入势公式；通过组件测试证明共享状态已传到三条路径。
- 保留 `RINTI(...,1)` 的动能语义，不在 `SETHAM` 中再叠加独立 NP 矩阵。
- 不增加新的二电子算符，不重写角向系数、CSF 生成、EOL 权重、冻结轨道或 SCF 停止条件。
- 不把 NP 初始化放进 `IMPROV`、`SOLVE` 或每次 `MATRIX` 调用中。
- 核对 `_mem` 的独立初始化及 `SETHAM` 文件，不能以非内存版本完成接入代替全部四个目标。

### 4.5 MPI 路径与错误协议

1. 所有 rank 先核对本地活动长度/容量；用集体操作归并失败状态，避免部分 rank 提前退出。
2. root 在 `permdir` 对应路径读取和验证配置，广播成功/失败状态及固定长度错误消息。
3. 正常时广播配置中的标量与字符串；使用显式字段广播，不把含 allocatable 字段的派生类型直接作为原始字节广播。
4. root 在最终径向网格上构造 NP 数组；再次广播构建状态，成功后广播 `N` 及活动数组。
5. 各 rank 从各自刚生成的纯 Coulomb 数组和同一份 NP 数组合成完整 `ZZ`。
6. 新配置错误路径在全部 rank 同步后统一返回失败；需要立即终止时调用 `MPI_Abort` 并给出非零码。不得仅 root 调用含 barrier 的 `stopmpi`，也不能以普通 `stop` 的零退出码掩盖失败。

继续使用原有 `YP=ZZ/NPROCS` 与 `COFPOTmpi` 求和归约。**不要除小存储的 `ZZ` 或 NP 广播数组**；`DACON` 与 `RINTI` 仍需要完整外势。仅新增 NP 特定适配器，不改写现有通用 MPI 启动/目录帮助函数。

### 4.6 输出与兼容性

保留 `rwfn.out`、`rmix.out` 和旧输入的二进制/交互格式。启用 NP 时由串行进程或 MPI root 写一份独立 `npol.summary`，内容限于：配置版本、通道开关、核标识、参数来源、内部单位、强度/截断参数、`N`、是否成功完成势合成。

默认关闭且文件不存在时不新增交互问题、不写 NP 专属文件、不改变已有提示文本。显式关闭的配置也不能改变数值数组。输出打开/写入失败必须返回状态，MPI 侧同步错误，不能形成假成功摘要。

## 5. 代码测试设计

### 5.1 测试边界

全部新增测试放在 `test/npol/`。测试输入仅为人工标量/数组、临时 namelist 文件和小型程序桩；调用函数或组件后立即退出，不启动完整 `rmcdhf` 作业，不依赖真实 `isodata`、`rwfn.inp`、`mcp.*` 或外部准备工具。

数学求值断言属于内核代码单元测试。MPI 测试只验证配置广播、数组一致性、归约和进程退出，不执行电子轨道优化。本文不注册或调用现有的原子计算回归、SLURM 脚本及生产数据目录。

### 5.2 测试矩阵

| 拟注册 CTest 名 | 测试方法 | 必须捕获的问题 |
| --- | --- | --- |
| `npol_code_input` | 临时文件覆盖缺失文件、关闭、单通道、双通道、未知字段/版本/单位、缺参、NaN/Inf、非法截断、错误核标识和 I/O 失败 | 解析失败被误当关闭、缺参被默认零吞掉、单位重复转换 |
| `npol_code_kernel` | 人工常量与半径，检查 $r=0$、$r=b$ 的闭式值、通道开关、显式零强度、叠加及有限性 | 分母幂数错误、E1/E2 串用、重复乘系数、关闭通道仍求值 |
| `npol_code_compose` | 固定 Coulomb 哨兵数组，合成两次、关闭/启用切换、重新 reset 后换网格 | 加号/减号错误、重复叠加、旧状态泄漏、把总势当基线 |
| `npol_code_bounds` | `N=1`、`N=NNNP`、非法长度、活动范围外哨兵；GNU 调试构建启用数组边界检查 | 越界、尾部被写坏、原点除零、非法网格未拒绝 |
| `npol_code_consumers_<variant>` | 分别编译四个变体的生产 `YPOT/DACON` 和公共 `RINTI` 到独立测试驱动，注入人工状态 | 只接入直接势、遗漏非对角项、`MODE=1` 污染、对共享势重复加项 |
| `npol_code_upstream_copy` | 剥离明确列举的 NP hook 后检查 492 个复制文件的 SHA-256 | 误用 orbopt 源码、混入其求解流程或丢失原版文件 |
| `npol_code_entrypoints` | 对四个入口、接口和主调用点作结构检查；结合生产 `NPOL_INIT` 的动态驱动 | 漏接 `_mem`、签名不同步、hook 在 `NUCPOT` 前或轨道读取后、MPI 路径依赖临时 cwd |
| `npol_code_mpi_broadcast_1/2` | 用真实 MPI 启动独立 NP 测试驱动，root 持有配置，worker 工作目录没有配置文件 | worker 意外 I/O、字段未广播、rank 间数组不同、错误的进程数缩放 |
| `npol_code_mpi_error_2` | root 解析失败以及非 root 长度错误两个子场景，父进程检查非零退出和有界终止 | 单 rank `STOP` 导致挂起、失败返回零码、先失败后仍广播数组 |
| `npol_code_output` | 临时目录中写回/读取元数据摘要；关闭时检查无输出 | 单位/通道标记与实际数组不符、worker 重复写文件、I/O 错误被忽略 |

`<variant>` 对应四个明确支持目标，测试驱动按独立可执行文件编译，避免不同应用的同名 Fortran 例程在同一个测试中冲突。MPI 缺失时不能伪装这些项目已通过；报告未构建/未运行。

### 5.3 消费者测试的独立预期

不要只测试新 helper，必须覆盖实际生产消费者。可在测试可执行文件中使用确定的 `DPBDT/QUAD` 桩，隔离外势路径；这些桩仅属于 `test/npol/`，不链接进正式程序。

- **YPOT：**令 `NYCOF=0`，串行返回完整 `ZZ`。MPI 变体返回 rank-local `ZZ/NPROCS`；小型 MPI harness 求和后应恢复一份完整 `ZZ`。
- **DACON：**令 `NDCOF=1`，固定系数 `q` 及人工 `PF/QF`，比较只改变 `ZZ` 时的输出差。预期为 $\Delta XP_i=-q\,\Delta ZZ_i\,Q_i/C$、$\Delta XQ_i=q\,\Delta ZZ_i\,P_i/C$。这直接检查 `ZZ/C` 路径及符号。
- **RINTI：**固定导数和带权求和桩，使预期值可手工计算。改变 `ZZ` 时，`MODE=0` 的变化应是 $-\sum_i w_i\,RPOR_i\,\Delta ZZ_i(P_{ia}P_{ib}+Q_{ia}Q_{ib})$，`MODE=1` 应不变。另放置仅非对角轨道对非零的夹具。
- **关闭兼容：**比较接入前后同一人工状态的所有消费者输出。对无需归约的复制/关闭路径检查精确相等；对浮点求和使用按机器 epsilon、操作数尺度和求和项数确定的容差，不借用物理误差阈值。
- **初始化组合：**动态调用正式 `NPOL_INIT`，确认成功才提交新状态、失败不留下部分总势。hook 结构检查只能补充动态测试，不能代替上述消费者覆盖。

内核预期值使用简单闭式值，例如 $V(0)=-K/b^p$、$V(b)=-K/(2b^p)$，其中 $p=4$ 或 $6$。测试中选用的人工数值必须标为软件夹具，不描述为任何同位素参数。

## 6. CMake 与测试运行入口

修改既有 [libmod CMake](../../src/lib/libmod/CMakeLists.txt)、[lib9290 CMake](../../src/lib/lib9290/CMakeLists.txt)、[mpi90 CMake](../../src/lib/mpi90/CMakeLists.txt) 和 [test CMake](../../test/CMakeLists.txt)。应用源列表仅登记它们实际需要的新本地包装文件；共享实现只编译一次。

- 使用现有 `setup_fortran_modules`、`target_link_libraries_Fortran` 管理 `.mod` 搜索路径。
- 给所有新增测试设置统一 `LABELS npol_code`。MPI 测试还设置 `PROCESSORS`、有限 `TIMEOUT`，并按 `TARGET mpiu90` 与 `MPIEXEC_EXECUTABLE` 条件注册。
- MPI 启动使用 CMake 提供的 `${MPIEXEC_EXECUTABLE}`、`${MPIEXEC_NUMPROC_FLAG}` 和 pre/post flags，不硬编码某一种 MPI 命令。
- MPI 失败测试由测试驱动包装器判断预期非零状态；超时本身必须判为失败，不能简单用 `WILL_FAIL` 将所有异常当成成功。
- 如需 Python 包装器，沿用 `../graspkit-tools/.venv/bin/python` 和现有 `Python3_EXECUTABLE` 选择机制；不创建新环境、不执行依赖同步。
- 不把需要 `GRASP_PREPARATION_BINDIR` 或实际原子数据的既有测试加到 `npol_code` 标签。

实施后使用新的构建目录，以下代码验证命令已于 2026-10-24 执行，结果见文末：

```bash
# 从 grasp-atomic-suite 仓库根目录运行。
cmake -S . -B build-npol-upstream-serial -DCMAKE_BUILD_TYPE=Debug \
  -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-npol-upstream-serial -j4
ctest --test-dir build-npol-upstream-serial -L '^npol_code$' --output-on-failure

cmake -S . -B build-npol-upstream-mpi -DCMAKE_BUILD_TYPE=Debug \
  -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-npol-upstream-mpi -j4
ctest --test-dir build-npol-upstream-mpi -L '^npol_code$' --output-on-failure
```

先确认 MPI 配置确实找到工具链、四个支持目标已构建、对应测试已注册。串行配置应只注册非 MPI 项；MPI 配置中的缺项必须报告。GNU 编译器下可对新测试目标添加 `-fcheck=all -fbacktrace`；保持既有公共编译选项和默认容量不变。这里不执行 `install`，也不运行不带过滤条件的全量 CTest。

## 7. 实施顺序与代码完成条件

0. **原版复制**：撤回旧 orbopt 接入，从本地 grasp 复制四个 rmcdhf 应用，记录源码摘要；保留原版求解流程和旧 orbopt 文件。
1. **输入和内核**：确定上述 v1 契约，新增类型、解析、校验、纯求值与合成例程；完成 input/kernel/compose/bounds 测试。
2. **串行接入**：接通两个串行入口和接口，保留三个生产消费者；完成串行消费者、默认关闭和输出测试。
3. **MPI 接入**：显式传递 `permdir/npol.inp`，加入状态同步和广播，接通两个 MPI 入口；完成独立 MPI 协议及错误退出测试。
4. **构建和文档收尾**：确认四个新原版目标及测试标签、旧 orbopt 无 NP 接入，记录支持矩阵和输入格式，更新文档索引；执行上述有范围的代码测试。

完成时逐项记录：

- [x] 四个应用来自 grasp 新复制的原版快照，已有 orbopt 文件恢复到本轮改动前版本。
- [x] 新输入的关闭、启用、错误行为有确定契约，旧交互输入无需更改。
- [x] 总势在初始化时只组合一次，重复合成与重新初始化的状态明确。
- [x] `YPOT`、`DACON`、`RINTI(...,0)` 的生产路径均有测试证据。
- [x] `RINTI(...,1)`、冻结轨道、EOL 权重、公共核模型/网格和未接入应用的默认行为不变。
- [x] MPI 的完整外势、局部归约、目录选择和集体错误处理均由独立测试驱动覆盖。
- [x] 支持目标全部构建；测试结果按“通过、失败、未运行”记录，不把缺失项计为通过。
- [x] `npol_code` 测试不依赖真实原子数据，不启动 SCF、不提交作业。
- [x] 现有二进制输出格式保持兼容，NP 摘要与配置一致。

以上完成条件只证明软件契约及接入正确性，不作物理精度结论。实施完成后再将现行接口说明写入 `docs/implemented/`，本计划保留完成状态与实现提交，避免将设计文字误读为已实现功能。


## 8. 实施结果（2026-10-24）

- 接入基线已更正为从 grasp 新复制的四个原版应用；旧 orbopt 的 16 个改动文件已恢复，原版 492 个复制文件有 SHA-256 核对。
- 已新增配置/状态、输入校验、内核、串行初始化及 MPI 适配模块；四个新目标均在 NUCPOT 后显式接入。
- 串行 Debug 构建完成，`npol_code` 9/9 通过；MPI Debug 构建完成，四个支持目标均生成，`npol_code` 14/14 通过。
- 额外选择既有纯代码回归：串行 7/7、MPI 配置 10/10 通过；没有运行全量 CTest 或原子数据回归。
- 真实 MPI 检查包含 1/2 ranks、两个生产消费者变体、完整字段/数组广播、临时 cwd，以及七种集体失败场景的非零退出。
- 关闭路径、生产 YPOT/DACON/RINTI、MODE=1 不变性、重置/重复合成、边界和摘要失败均有组件测试；冻结轨道、EOL、SCF 判据和旧二进制输出的实现未修改，未把此项源码不变性扩大为真实 SCF 回归结论。
- 输入具体上限为 65536 字节，字符串上限为 units 16、model_id 64、parameter_source 1024 字符；拒绝额外 namelist 内容。
- 在现有 `-fno-automatic` 下采用固定解析缓冲区，避免重复初始化遗留的 allocatable 状态；摘要增加写后内容校验以捕获运行库可能漏报的写错误。
- 公共模块/例程仍只由原有库提供；只有 MPI 适配模块依赖 MPI。Python 环境与共享容量未修改；本次提交仅包含源码、测试和文档，未安装或推送。

实现说明：[rmcdhf_nuclear_polarization.md](../implemented/rmcdhf_nuclear_polarization.md)。
测试入口与边界：[test/npol/README.md](../../test/npol/README.md)。
上述勾选仅表示本文规定的软件验证完成；核素参数、SCF 响应及物理精度验证仍未执行。
