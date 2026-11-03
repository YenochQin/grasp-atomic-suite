# jj2lsj_omp：普通 CSF 的 OpenMP jj–LS 转换

状态：2026-11-02 已实施；完成合成回归和 Gd I 真实输入兼容性验证，已有单次性能观测；不代表原子物理精度或收敛验收。
suite 修改前基线：`a4661441492e38beb99293182cb9257326d2dcef`。

## 来源与适用范围

普通 GRASP `jj2lsj90` 源码来自工作区的 `grasp` 快照，文件头最新标注 January 2025；
没有可确认的独立上游提交号，因此用[来源文件 SHA-256 清单](../reference/jj2lsj/source-manifest.json)
锁定导入前的内容。保留原作者与 MIT 许可；此清单不是集成后源码的哈希。

新程序是 `jj2lsj_omp`，读取 `isodata`、`name.c`、`name.m` 或 `name.cm`。
交互问题、阈值默认值、`.lsj.lbl/.lsj.c/.lsj.j`、唯一标签和 `.lsj.T` 格式延续来源程序。
已有 `jj2lsj_csfg` 仍是独立 CSFG 应用，本次没有修改它的转换算法。

## 构建与运行

从 suite 仓库根目录执行：

```sh
cmake -S . -B build -DGRASP_JJ2LSJ_OPENMP=ON
cmake --build build --target jj2lsj_omp test.jj2lsj_coefficients --parallel 4
ctest --test-dir build -R '^jj2lsj_' --output-on-failure
cmake --install build --component jj2lsj_omp
```

新目标默认构建，不依赖 MPI、CSFG 或 `GRASP_BUILD_SERIAL_APPS`。
`GRASP_JJ2LSJ_OPENMP` 默认为 `ON`；找不到 OpenMP Fortran 时明确警告并构建串行版本。
显式 `-DGRASP_JJ2LSJ_OPENMP=OFF` 可关闭 OpenMP。
已验证 Linux / GNU Fortran 16.2.1；其他编译器尚未验证。

在计算输入目录中运行已安装的程序，例如：

```sh
OMP_NUM_THREADS=8 /path/to/grasp-atomic-suite/bin/jj2lsj_omp < jj2lsj.inp > jj2lsj.log
```

`/path/to/grasp-atomic-suite` 是需要替换的安装位置，不是新的输入格式。
`OMP_NUM_THREADS=1` 使用保留的串行转换路径；也可设置 `JJ2LSJ_SERIAL=1`
在不改变线程环境的情况下强制走串行路径，以便同一构建做对照。
其他 `JJ2LSJ_SERIAL` 值按未启用处理。
并行路径在每个转换 block 开始时报告最大线程数与每批 LS 行数。
不需要用 `mpirun` 启动此程序。

## 并行与数值契约

实现位于 [jj2lsj_code.f90](../../src/appl/jj2lsj90/jj2lsj_code.f90)：

- `asf2ls_serial` 保留原来的变换循环，作为单线程及变换文件路径。
- `transform_ls_row` 完整计算一行 LS-CSF 的所有所选能级系数。每行保持原来的
  jj-CSF 求和次序、`EPSNEW` 和角动量计算，不引入占据预筛或浮点归约。
- `asf2ls` 用 OpenMP `parallel do` 动态分配批内各行。批次按原 LS 顺序提交，
  按原比较式累积权重并检查 `MINCOMP`；首个满足条件的行之后的预计算结果被丢弃，
  原有零尾部不写入。不同批次之间不提前提交或跨序执行截断判定。
- 每批最多 256 行，通常为最大线程数的四倍，并将系数缓冲限制在约 8 MiB。
  所选能级、累计权重、标志数组和线程缓存另有小额内存开销。
- `.lsj.T` 写入/回读（`ioutT=1/2`）始终回退串行，保持记录及索引顺序。
  因而这些模式目前不会获得 LS 行计算的 OpenMP 加速。
- 输入、`setLS` 基生成、排序、标签输出和 `uniquelsj` 都仍为串行。
  同一 block 的全部所选能级共同决定 `ithresh`，不按能级拆成独立转换。

为了隔离历史静态变量，[应用 CMake](../../src/appl/jj2lsj90/CMakeLists.txt) 仅在自身目录
去掉 `-fno-automatic`，GNU 构建加 `-frecursive`。串行控制例程用显式 `SAVE`
保留原工作区行为；变换函数和每行工作区可重入。

`jj2lsj_thread_kernels` 是一个 **object target**，从现有共享源码重编并发调用闭包：
`nine/nine0/sixj/sixj1/2/3/35/4/5/dracah/ittk/ixjtik/iq/jqs/itjpo/ispar`。
这些对象直接链接到调用程序，先于历史静态库解析同名符号。
没有复制 `libmod`、`lib9290` 或它们的全局状态；其他 suite 程序继续使用原库。
转换库的接口会将这些对象传递给测试驱动，避免测试误链接旧的静态内核。

共享 [sixj.f90](../../src/lib/librang90/sixj.f90) 只加入 `!$OMP THREADPRIVATE(RACA)`：
OpenMP 专用对象中的惰性缓存由各线程独占，普通 suite 构建将这条指令视作注释。
`FACTT` 在启动线程前初始化的阶乘表，以及输入数组、常量表，只读共享。

导入时同时处理三处局部生命周期/边界问题：释放 LS-CSF 的 `w` 数组；显式初始化
唯一标签循环的 `New_J`；在 [packLS.f90](../../src/appl/jj2lsj90/packLS.f90) 中使用有边界保护的
前一个占据数，避免 Fortran 非短路逻辑表达式访问 `Q(0)`。
后者由 `-fcheck=all` 检查发现；修复后的参考输出保持不变。

## 已完成验证

[回归说明](../../test/jj2lsj/README.md) 和 [测试脚本](../../test/jj2lsj/check.py)
记录合成输入与预期结果的来源。测试不需要 SCF、RCI、外部 Python 包或集群作业。

1. Release OpenMP 构建及 `jj2lsj_parallel_regression`：48 组系数/输出比较，
   覆盖 1/2/4/4 线程配置（4 线程冷启动重复），逐字节比较完整带符号双精度系数及输出文件。
2. Debug `-g -fcheck=all -fbacktrace` 构建；同一回归，包括 OpenMP 关闭后的串行构建。
3. 完整变换的独立正交归一检查；默认截断、自定义截断、跨批次、非整批尾部、
   跳过 block、非顺序选择能级、唯一标签、`.lsj.T` 写入和回读。
4. 以原始 GRASP 可执行文件对照三个体系 × 四种选择/截断设置，
   `.lsj.lbl/.lsj.c/.lsj.j` 输出一致；默认标签参考保存在测试目录。
5. suite 原有 `lib9290_quad`、`shared_grid_capacity`、`csfg_capacity`、
   `csfg_file_formats`、`csfg_alternate_capacity` 五项回归通过。

另对仓库已有 `data/test.c`、`data/test.m` 做默认模式回放，原始 GRASP、新串行路径、
4 线程结果的标签文件一致。此样例三条路径均报告相同的 `IEEE_DIVIDE_BY_ZERO`
浮点状态提示；它是对照中已有的现象，未作为本次并行实现的物理验证结果。

同日以提交 `3630f42` 测试 Gd I 的六套真实输入，18 个输出文件与历史结果一致。
ZF-AS1、RCI-AS4 共 994,894 个默认筛选后的 LS 系数在 4/8 线程下与串行逐位一致，
包含 AS1 的实际提前停止分支。单次转换核加速约为 4 线程 3.8×、8 线程 5.5–5.9×；
完整设置、端到端时间及既有提示见[真实输入验证报告](../reviews/jj2lsj_gd_i_validation_2026-11-02.md)。

## 限制与后续

当前只有单次性能观测，没有重复统计的性能基准或真实原子精度/收敛结论。合成回归及上述真实输入
验证仅覆盖所测路径，不能推广到所有组态、编译器或线程配置；合成测试使用人为构造的混合向量。

总体速度仍可能受串行 `setLS` 去重/递归生成、排序、I/O 和唯一标签限制。
输入 `EVEC` 与 LS 结果仍沿用原来的全局分配方式，唯一标签也保留原固定容量。
未实施 MPI、大规模数据布局压缩、占据数预筛或重新排序。
下一步应在固定环境中重复计时，并进一步拆分输入和排序的开销，再决定是否扩展优化。
按 ASF 并行后处理、唯一标签协调和跨 block 异步流水线的待执行讨论见
[后处理并行与异步方案](../plans/jj2lsj_postprocessing_parallelization_2026-11-03.md)。
