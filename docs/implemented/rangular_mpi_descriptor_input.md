# MPI rangular 描述符输入与性能对照

日期：2026-11-11。实现基于 `c47f4d42b1f60bd77eb04dd3a243f545d124f481` 的工作树增量。
普通 `rangular_mpi` 已支持与串行/OpenMP 相同的 GCSDESC1 二进制；默认仍使用文本 `rcsf.inp`。
本次未改变角系数内核、MPI 列分工、系数排序及 MCP 分片格式。CSFG、RMCDHF 和 RCI 的输入入口不在本次范围内。

## 输入适配

[rangular_mpi_input.f90](../../src/appl/rangular_mpi/rangular_mpi_input.f90) 负责选项、路径及广播。
描述符路径由 rank 0 调用现有 [csf_descriptor_input](../../src/lib/lib9290/csf_descriptor_input.f90) 解码，
安装头部后广播 block/轨道元数据，每个 block 装配原 IQA/JQSA/JCUPA 数组并广播给其它 rank。
其它 rank 无需打开二进制文件，也不持有描述符原文件的副本。
byte 数组通过 `MPI_INTEGER1` 发送，按完整 CSF 列拆分消息，限制单次 MPI 默认整数 count，避免广播长度乘法溢出。
这一保护没有扩大其它科学内核的整数容量。

文本路径继续调用原 `CSLHMPI/LODCSLMPI`。两条路径沿用相同的相互作用设置和后续计算。
新变量使用实际含义名称，旧数组名仅在 `USE ... ONLY: occupations => IQA` 等连接处保留。
文件协议、支持的物理态和预处理方法见[公共描述符说明](rangular_descriptor_input.md)。

| 选项 | 行为 |
| --- | --- |
| 无输入选项 | rank 0 从 MPI serial I/O 目录读取 `rcsf.inp`。 |
| `--csf-text FILE` | 指定原文本文件；与描述符选项互斥。 |
| `--csf-descriptors FILE` | 指定预处理二进制；无需另外提供 `rcsf.inp`。 |
| `--timings timings.csv` | 各 rank 独立生成 `timings.rank000.csv`、`timings.rank001.csv` 等；默认关闭。 |

相对输入路径以 rank 0 的 serial I/O 目录为基准，通常为启动目录，也可以由旧 `disks` 文件配置。
相对计时路径以各 rank 的工作目录为基准；绝对计时路径也追加 rank 后缀，避免进程间覆盖。
计时文件必须不存在。输入/计时文件名拒绝 `mcp*`、`tmp*`、`rangular.log`、`genmcp.dbg` 等保留名称。
绝对输入路径使用可分配字符串；MPI 启动工作目录及 `MPI_TMP` 仍受旧 CPATH 的短路径约束。

## 构建与运行

从 suite 根目录使用 fish：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
cmake -S . -B build-mpi-descriptor \
    -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
    -DCMAKE_BUILD_TYPE=Release -DGRASP_FORTRAN_PROFILE=f2018-release \
    -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
    -DGRASP_BUILD_CSFG=OFF -DGRASP_RANGULAR_OPENMP=OFF -DGRASP_JJ2LSJ_OPENMP=OFF
cmake --build build-mpi-descriptor --target rangular_mpi csf_descriptor_verify -j4
```

在独立、路径较短的计算工作目录运行，例如：

```fish
set -gx MPI_TMP (pwd)
mpirun --bind-to none -np 8 /path/to/suite/build-mpi-descriptor/bin/rangular_mpi \
    --csf-descriptors /path/to/csfs.bin --timings timings.csv
```

完整/一阶相互作用仍由原交互流程选择。MPI 原有一阶输入要求每个 block 满足 `1 < ICCUT < NCFBLK`；
本次没有改变此约束。输出仍位于 rank 工作目录，文件名如 `mcp000.30`、`mcp001.31`。

## 可重复的文本/描述符对照

[compare_rangular_mpi_inputs.py](../../scripts/compare_rangular_mpi_inputs.py) 使用同一个 MPI 可执行文件：

1. 固定输入副本并保存文件/程序 SHA-256；先用原生验证器逐项核对描述符与文本数组。
2. 每轮交替运行文本和描述符，固定数学库为单线程，完整相互作用由 stdin `y` 指定。
3. 保存每个 rank 的计时及完整 MCP，逐文件比较整个输出集合与第一轮文本基准。
4. 统计进程外完整墙钟和 rank 0 的加载/广播总时间。输入准备、数组验证和 MCP 哈希比较均不计入计时。

从 suite 根目录运行，`--output` 必须是新目录：

```fish
set result_directory temp/mpi-descriptor-(date +%Y%m%d-%H%M%S)
../graspkit-tools/.venv/bin/python scripts/compare_rangular_mpi_inputs.py \
    --executable build-mpi-descriptor/bin/rangular_mpi \
    --verifier build-mpi-descriptor/bin/csf_descriptor_verify \
    --csf /path/to/rcsf.inp --descriptors /path/to/csfs.bin \
    --ranks 8 --repeat 3 \
    --mpiexec-arg=--bind-to --mpiexec-arg=none \
    --output "$result_directory"
```

`--mpiexec` 和 `--numproc-flag` 可指定 launcher；额外 launcher 参数用可重复的 `--mpiexec-arg=VALUE`。
`--scratch-root` 默认使用系统临时目录，内部生成短的 rank 工作路径，避免较长的结果目录触发旧 CPATH 限制。
所有计时完成后才将运行目录移入结果目录，避免归档输出扰动下一轮；因此 scratch 必须能容纳全部重复的输出。
在服务器测量磁盘性能时，应显式把 scratch 放在实际计算所用文件系统；默认临时目录可能是 tmpfs。
不要同时编译、启动其它计算或执行大规模数值比较。

`summary.json` 的 `status=passed` 表示数组验证、每个 rank 的完整计时和全部 MCP 字节比较通过，
不表示性能一定提升。`statistics.by_mode` 保存原始时间及中位数/范围，`statistics.ratios` 保存相对变化；
正的 `median_reduction_fraction` 表示描述符较快，负值表示较慢。
`array-verification.log`、各轮日志、CSV 和输出都保留。原始 descriptor 生成/预处理成本不包含在运行时间内。

## 计时边界与验证

MPI 复用 [rangular_profile](../../src/appl/rangular/rangular_profile.f90) 的 CSV 格式，不新增逐 CSF 对计时。
每个 rank 的计时相互独立，没有为了记录时间添加全局屏障或归约。

| 阶段 | 内容 |
| --- | --- |
| `input_header_and_broadcast` | 文本扫描或二进制头部加载及广播。 |
| `csf_allocation` | 本 rank 的 block 运行数组分配。 |
| `csf_load_and_broadcast` | rank 0 加载/语义装配及数组广播；其它 rank 还可能包含等待前一阶段完成的时间。 |
| `mcp_setup_including_interaction_input` | MCP 设置及可能的一阶交互等待。 |
| `pair_generation_and_buffering` | 本 rank 的角系数计算及临时记录写入。 |
| `sparse_structure_and_output`、`sort_and_output` | 本 rank 的结构、排序与输出；排序按通道记录。 |
| `program_before_finalize_including_interaction_input` | `startmpi2` 返回后至 `stopmpi2` 调用前；不含 MPI 启动/结束。 |

运行器的完整墙钟包含 MPI 启动与结束。rank 0 的输入指标包含广播和同步影响，不能称为纯文件读取时间。
不同 rank 的阶段会重叠，各阶段最大值也可能来自不同 rank，不能直接相加代替整程耗时。

[六项 MPI 回归](../../test/rangular/test_mpi_descriptors.py) 在 Release 和 F2018 Check 下通过：
完整相互作用的 1/2/4 ranks、空任务 rank、一阶限制、相对/绝对路径、二进制运行不依赖文本文件、
无效文件/物理态的非零退出，以及计时文件防覆盖。
跨 rank 数使用 MCP 解析器零容差比较，同 rank 数的文本/描述符文件逐字节一致。

```fish
ctest --test-dir build-mpi-descriptor -R '^rangular_mpi_descriptors$' --output-on-failure
```

新适配器及共享计时模块在本次 MPI Check 编译中无警告。三个原有 labeled DO 警告仅移动行号；
警告预算按 MPI 目标增加两个编译单元，没有提高诊断数量预算。其它矩阵配置作结构性更新，没有重跑全部矩阵。
实际输入的性能及完整输出对照见 [2026-11-11 MPI 描述符测量](../reviews/rangular_mpi_descriptor_benchmark_2026-11-11.md)。
