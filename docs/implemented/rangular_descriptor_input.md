# rangular 描述符输入、数组对照及阶段计时

日期：2026-11-10。实现基于 suite `e229657d8a32b851da9c130928545bcb3ac209f1` 的工作树增量；
描述符契约核对 rCSFs `6e4fcff1439fb1b56d4891908c5bd58722ddf64a`，实际兼容性检查使用已安装的原生 rCSFs 3.0.0。

状态：普通串行 `rangular` 已接入；原文本入口保留。预处理、原文本/新输入数组对照、小型 MCP 回归和阶段计时已实施。
后续增加了独立的 [rangular_omp](rangular_openmp.md)；本文仍描述两种输入的公共契约。
2026-11-11 已增加普通 [rangular_mpi 描述符入口](rangular_mpi_descriptor_input.md)，由 rank 0 解码并广播。
CSFG/RMCDHF/RCI 的输入入口未改造。
这些结果验证软件和输入表示，不建立生产规模性能、活性空间收敛或物理精度。

后续基于 `c5f36d2` 的 [2026-11-11 本机性能对照](../reviews/rangular_mpi_openmp_benchmark_2026-11-11.md)
已完成 14,551-CSF 样本的串行/OpenMP 文本与描述符测量：串行输入阶段中位数约快 26.3 倍，
整程收益较小；该报告单独列出转换成本和样本限制。MPI 的后续实现与测量见上一段链接。

## 1. 代码归属与运行路径

| 文件 | 职责 |
| --- | --- |
| [prepare_csf_descriptors.py](../../scripts/prepare_csf_descriptors.py) | 分批读取 Parquet 与 header TOML，校验版本/类型/绑定/分块/基本范围，输出版本化二进制。 |
| [csf_descriptor_input.f90](../../src/lib/lib9290/csf_descriptor_input.f90) | 解码交换格式；复用 `terms_C` 校验物理态；补全 J/K、组织数组；提供头部及逐 block 加载接口。 |
| [rangular_input.f90](../../src/appl/rangular/rangular_input.f90) | 选择原 `CSLH/LODCSH2` 或描述符路径；连接现有程序状态。 |
| [csf_descriptor_verify.f90](../../src/appl/rangular/csf_descriptor_verify.f90) | 独立调用两条加载路径，逐项比较元数据及全部三个 byte 数组。 |
| [verify_csf_descriptors.py](../../scripts/verify_csf_descriptors.py) | 运行原生验证器并要求成功标记，防止旧文本 reader 的 `STOP` 返回零时被误判通过。 |
| [rangular_profile.f90](../../src/appl/rangular/rangular_profile.f90) | 可选的 SYSTEM_CLOCK 墙钟计时，按阶段写 CSV。 |

新增变量采用 `principal_numbers`、`relativistic_kappas`、`subshell_quantum_numbers`、
`coupled_two_j`、`block_index` 等含义明确的名称。原模块名通过 `USE ... ONLY: 新名称 => 原名称`
在边界连接；没有整体重命名旧角系数计算代码。

数据路径为：

```text
描述符 Parquet + 匹配的 header TOML
  → suite Python 预处理
  → GCSDESC1 二进制
  → Fortran 头部与 block loader
  → 原 IQA / JQSA / JCUPA 和原串行角系数流程
```

预处理需要 Python >=3.11、NumPy 和 PyArrow；本轮使用共享环境 Python 3.15.2、PyArrow 26.0.0。
suite 的编译和 rangular 运行不依赖 Python、PyArrow 或 Rust。维护测试中的 Parquet 夹具自包含，不调用 rCSFs。
严格串行验收现在要求该集成测试，因此测试用 Python 环境须提供 NumPy/PyArrow；缺依赖会明确失败，不静默跳过。

## 2. 构建和使用

需要直接运行原版和新版并比较结果时，使用[版本对照运行器](rangular_version_comparison.md)。
它组合预处理、数组验证、独立目录运行和 MCP 比较，下面保留各工具的独立用法。

以下构建命令从 suite 根目录执行，使用独立目录，不安装或覆盖现有程序：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
cmake -S . -B build-descriptor-check \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
  -DCMAKE_BUILD_TYPE=Debug -DGRASP_FORTRAN_PROFILE=f2018-check \
  -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON \
  -DGRASP_BUILD_CSFG=OFF -DGRASP_JJ2LSJ_OPENMP=OFF
cmake --build build-descriptor-check --target rangular csf_descriptor_verify -j4
```

从已有、与 header 绑定的完整 V2 描述符生成输入；路径只是使用示例：

```fish
../graspkit-tools/.venv/bin/python scripts/prepare_csf_descriptors.py \
  /path/to/descriptors.parquet /path/to/csfs.bin \
  --header /path/to/source_header.toml

../graspkit-tools/.venv/bin/python scripts/verify_csf_descriptors.py \
  /path/to/csfs.bin /path/to/original.c \
  --verifier build-descriptor-check/bin/csf_descriptor_verify
```

预处理默认每批 8192 行，可用 `--batch-rows` 调整；它输出 JSON，包含行数、block 数、文件大小、
header SHA-256 和预处理墙钟时间。输出文件必须尚不存在；转换中捕获的异常会删除本次不完整输出。
命令行输入/格式错误退出 2。验证 wrapper 通过退出 0，失败退出 1，并显示首个差异的位置。
不要只根据原生验证器的退出码判断旧文本输入是否通过，必须同时要求其最终成功标记；wrapper 已执行此检查。

在独立计算工作目录中运行编译后的程序，使用其绝对路径。例如：

```fish
/path/to/suite/build-descriptor-check/bin/rangular \
  --csf-descriptors /path/to/csfs.bin --timings timings.csv
```

无参数仍读取当前目录的 `rcsf.inp`；`--csf-text FILE` 可指定其他文本路径。
`--csf-text` 与 `--csf-descriptors` 互斥。二进制路径不需要工作目录中存在 `rcsf.inp`；
完整/一阶相互作用及各 block 的 `ICCUT` 仍按原交互流程输入。
计时文件必须尚不存在。输入和计时路径拒绝 `mcp.*`、`tmp.*`、`rangular.log`、`genmcp.dbg`
这些程序保留输出名称，避免打开同一文件进行相冲突的读写。

输出仍为串行 `mcp.*`。现有多 rank 求解器的 MCP 分片协议没有改变。

## 3. 支持与校验边界

预处理要求 Parquet **文件 footer** 中的 V2 元数据，而非只读取 `ARROW:schema` 的附加 metadata：
Rust writer 可以只在 footer 写入格式契约。必须为 physical seniority、非归一化、非空值 Int8 列，
字段顺序为每个 peel 壳层 `(n,2J,v,2K)`，尾部为 `(total_two_j,parity)`。
`source_header_sha256` 必须存在且与给定 TOML 字节完全相同；header-less 中间产物不作为可交付输入。
原始 rCSFs descriptor 文件不被修改。

保留 core/peel 次序、原行序、全部显式 block lengths，包括相邻同对称性 block。
不按 J/宇称重新排序，不改变 `ICCUT` 参考区间。检查电子数、block 内对称性、轨道重复、
占据上界、角动量范围、空壳的缺省通道、三角条件、耦合完整性和宇称。

Fortran 的数组语义：

- `IQA` 保存占据；固定 core 按满占据补入。
- `JQSA` 三槽为 `(physical seniority, w, 2J+1)`。仅在支持集合中将 `w` 置零。
- `JCUPA` 按非空非满的开放壳层链保存 `2K+1`；部分占据的 J=0 壳层仍保留对应耦合槽。
- `JCUPA(NNNW,:)` 保存 `(total_two_j+1)*parity`。两个 loader 都初始化全部数组，验证器比较包括零 padding 在内的全部值。
- 缺省 J 仅在给定占据/seniority 唯一决定物理态时补全；缺省 K 仅在零角动量传播、首个非零壳层或最终 J 已确定时补全。

普通文本 reader 之外的态不在本次扩展范围内：开放壳层单电子 `2j<=17`，
`2j>9` 仅允许电子/空穴数不超过 2；其余依现有 `terms_C` 表及原 reader 能区分的态检查。
`j=7/2,q=4` 的 seniority 分辨已验证。需要额外 `w`，或其他无法由普通 reader 区分的同 J 态，明确拒绝。
错误会报告 block 和 CSF 序号。`2J/2K<=126` 在窄化前检查；没有扩大原科学内核的整数容量。

当前 rangular 最多 50 个 block，轨道数受配置的 `NNNW` 限制；总 CSF 数和每 block 数须适合 default INTEGER。
交换格式的 64 位计数不代表求解器已经支持任意大基组。
本次没有实现全空间重复 CSF 检测，也没有修改旧 reader 的查重逻辑；生产输入仍须由合法生成/筛选流程提供。
格式检查不等于证明任意外部文件包含正确的科学基组。

## 4. 二进制 v1 协议

全部多字节字段为 **little-endian**；Fortran 按 byte 显式解码头部，不使用编译器特有 `CONVERT=`。
文件无 sequential-unformatted 记录标记；描述符量子数使用有符号 int8。

| 字节偏移 | 字段 | 类型/长度 |
| --- | --- | --- |
| 0 | magic `GCSDESC1` | 8 ASCII bytes |
| 8 | 版本，固定 1 | uint32 |
| 12 | 端序标记 `0x01020304` | uint32 |
| 16 / 20 | core / peel 轨道数 | 各 uint32 |
| 24 / 28 | block 数 / 每行通道数 `4*peel+2` | 各 uint32 |
| 32 | 总 CSF 数 | uint64 |
| 40 | 总电子数 | uint32 |
| 44 | 原 header SHA-256 小写十六进制 | 64 ASCII bytes |
| 108 | core 在前、peel 在后的轨道表 | 每项 `int32 principal_number, int32 kappa` |
| 轨道表之后 | 原顺序 block 表 | 每项 `uint64 row_count, int32 total_two_j, int32 parity` |
| block 表之后 | 原顺序描述符行 | 每行 `4*peel+2` 个 int8，行内按通道顺序连续 |

文件总长必须等于 `108 + 8*orbital_count + 16*block_count + row_width*total_csf_count`；
拒绝截断、尾随数据、不支持版本、非法端序标记、计数/尺寸不一致和窄化溢出。
Fortran 每批最多读 8192 行到临时数组，随后填充调用方当前 block 的三个原运行数组。
没有同时保留整份解码文件与运行数组。

header digest 标识预处理使用的 header；Fortran 不再打开 TOML，也不将该摘要当作 payload 校验和。
应保留输入文件及其生成记录，回归时仍用原文本身份与数组对照建立对应关系。

## 5. 计时语义

`--timings NEW_FILE.csv` 默认关闭，启用时写入：
`stage,block,channel,wall_seconds,records,storage`。
使用 int64 `SYSTEM_CLOCK` 计算墙钟，阶段边界写入，不在每个 CSF 对计时或写日志。
`block=0` 表示全局阶段；`channel=31` 为 T，`32+k` 为 V(k)，`-1` 表示不适用。

| 阶段 | 包含的工作 |
| --- | --- |
| `input_header` | 文本头部/block 扫描或二进制头部加载，storage 标出输入路径。 |
| `csf_allocation` / `csf_load` | 当前 block 数组分配；输入加载、语义验证和数组装配。 |
| `coefficient_buffer_allocation` | block 累积区和当前对缓冲的初始分配。 |
| `pair_generation_and_buffering` | CSF 对遍历、角系数计算、缓冲扩容和临时写入；records 为 T/V 贡献总数。 |
| `csf_and_pair_buffer_release` | 输入数组和当前对缓冲释放。 |
| `sparse_structure_and_output` | SETSDA 和稀疏结构输出；records 为稀疏位置数。 |
| `sort_and_output` | 每通道的 SORTMEM 或 SORT 及文件输出；storage 标出实际 memory/disk 分支，records 为该通道贡献数。 |
| `setup_including_interaction_input` | 建 MCP 头、设置相互作用参数及表初始化，可能含人工等待。 |
| `program_total_including_interaction_input` | 程序整体区间，包括人工输入等待和子阶段；不能与子阶段相加。 |

当前排序与输出按既有过程合并计时，成对计算与缓冲/临时 I/O 也没有在内循环拆开；
这些字段名明确反映边界，不将它们报告为纯计算或纯 I/O 时间。
仍有启动/日志/文件关闭等未单列开销，因此阶段之和不是严格的总时间分解。
计时开关下输出一致性已检查；本轮未建立生产规模计时开销或加速比。

## 6. 验证记录与复现

维护回归入口：

```fish
../graspkit-tools/.venv/bin/python test/test_csf_descriptors.py \
  --rangular build-descriptor-check/bin/rangular \
  --verifier build-descriptor-check/bin/csf_descriptor_verify -v
ctest --test-dir build-descriptor-check -R '^csf_descriptor_input$' --output-on-failure
```

[14 项测试](../../test/test_csf_descriptors.py)覆盖压缩/字典 Parquet、footer-only metadata、分批跨 block、
header 绑定与非法通道/格式拒绝、空/满/开放 J=0 壳层、缺省单壳 J、显式零耦合、
`j=7/2,q=4` seniority、相邻同对称性 block、`2J=126` 边界、额外 w 拒绝、差异定位与旧 STOP 防误判。
小型 rangular 回归有非零 T/V 贡献，覆盖完整相互作用和一阶限制。
文本/二进制输入 × 计时关闭/开启的八次运行，在同一构建内 MCP 文件逐字节一致；
独立 [MCP 解析器](mcp_compare.md) 同时检查结构与贡献。二进制运行时工作目录确实没有 `rcsf.inp`。

使用原生 rCSFs 3.0.0 对仓库现有 `test/jj2lsj/{even,odd,f3}.c` 实际转换并生成描述符后：

| 输入 | CSF 数 | block 数 | 与文本路径完全相同的数组值数 |
| --- | ---: | ---: | ---: |
| even.c | 82 | 7 | 52,070 |
| odd.c | 60 | 7 | 38,100 |
| f3.c | 29 | 5 | 18,415 |
| 合计 | 171 | 19 | 108,585 |

本地完整记录在忽略的 `build-descriptor-check/audit/rcsfs-native/`：各 case 保存原生生产器日志、
header/Parquet、交换文件和验证日志。此附加核查不成为维护测试对外部 rCSFs 的依赖。
本轮汇总、源码 SHA-256 及 CTest 日志保存在 `build-descriptor-check/audit/descriptor-input/`，
汇总文件为 `validation.json`；这些运行产物不纳入版本控制。

GNU Fortran 15.3.0：Check 使用 O0、运行时检查和浮点异常捕获；Release 使用 O2、严格接口和保守浮点选项。
串行/无 CSFG/无 OpenMP 的 Check 构建及 42 项 CTest 通过。
Release 全开构建通过；95 项测试中初次 39 项 MPI 测试被沙箱 socket 限制阻断，随后在允许本机通信的环境全部通过。
本次新模块没有编译器诊断；旧 rangular 的六条既有诊断仅因插入代码移动行号。

警告预算更新仅涉及四个新文件的编译覆盖和上述六个既有诊断的位置，数量未放宽。
实际测量 `check-s1-m0-c0-o0` 与 `release-s1-m1-c1-o1`；其他组合只按确定的 CMake 文件/目标增量更新预算，
**没有重跑全部 32 配置矩阵**，也没有消除既有 correctness-review 告警。
详见[警告预算](../reference/f2018/gfortran15-warning-baseline.json)中的 `descriptor_input_budget_update`。

未验证：百万级数据的性能/峰值内存、强制内存不足磁盘回退的计时路径、全量态空间、
CSFG 描述符入口，以及基于新输入的完整 SCF/RCI 科学结果。后续 OpenMP/MPI 的独立状态见文首链接。
