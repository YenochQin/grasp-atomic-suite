# RHFS_MPI `HFSGG_MPI` 本次优化说明

> 文档状态：历史第一阶段报告、部分描述已被后续实现覆盖。
> 核查日期：2026-10-05；源码基准：`8b298a6`。
> 外层筛选仍在，但文中计数器/细计时已在 3d067e7 删除，占据数筛选现用内联 IQA 扫描；后续已有零贡献筛选。正文耗时和瓶颈结论仅对应当时样例，不代表当前版本。
> 全部文档状态见[分类索引](../README.md)。

## 1. 优化范围

本次优化只涉及 `src/appl/rhfs90/hfsgg_mpi.f90` 中 `HFSGG_MPI` 的分布式主计算阶段，不改动物理公式、矩阵元定义和最终输出格式。

相关提交：

- `f008da6` `Add HFSGG MPI profiling detail`
- `cf197a2` `Pre-filter HFSGG MPI one-particle calls`

## 2. 优化前的主要问题

优化前，`HFSGG_MPI` 的总瓶颈不在 MPI 归约，而在内核循环中对 `ONEPARTICLEJJ` 的大量调用。

一次代表性运行结果显示：

- `HFSGG_MPI kernel wall time (max rank): 2048.128 s`
- `HFSGG_MPI reduce wall time (max rank): 15.238 s`

这说明：

- 通信只占很小比例，MPI reduce 不是主瓶颈
- 主时间消耗在每个 rank 的本地计算

进一步加 profiling 后发现，优化前每个 rank 大约有：

- `active ~= 1.51e9`
- `ia0 ~= 1.11e9`
- `nzacc ~= 2.81e7`

也就是说，进入 `ONEPARTICLEJJ` 的工作项里，有很大一部分最后直接以 `IA = 0` 返回，没有形成任何有效壳层贡献；真正进入最终累加的项只占很小比例。

## 3. 问题定位过程

为定位瓶颈，先在 `HFSGG_MPI` 中增加了轻量统计和阶段计时，重点观察：

- `active / onepcalls`
- `nzacc / zskip`
- `ia0 / diag / offd`
- `ptrig / ppar / pocc`

profiling 结果表明：

- `ptrig = 0`
- `ppar = 0`
- `pocc ~= 1.11e9`
- `ia0` 在优化后变成 `0`

这说明此前绝大多数空调用并不是由三角条件或宇称筛掉，而是由 `ONEPARTICLEJJ` 内部的占据数检查 `ICHKQ1` 拦截。也就是说，真正有效的优化方向不是改 MPI，也不是改调度，而是把这类廉价筛选前移到 `ONEPARTICLEJJ` 调用之前。

## 4. 本次实际优化点

### 4.1 将廉价筛选前移到 `HFSGG_MPI` 外层

在 `HFSGG_MPI` 中引入：

- `USE itrig_I`
- `USE ichkq1_I`

对应代码位置：

- [hfsgg_mpi.f90](../../src/appl/rhfs90/hfsgg_mpi.f90)
- [hfsgg_mpi.f90](../../src/appl/rhfs90/hfsgg_mpi.f90)

在调用 `ONEPARTICLEJJ` 之前，先做三类前置判断：

1. `ITRIG(ITJPOC, ITJPOR, 2*KT + 1)`
2. parity 检查 `ISPAR(IC)*ISPAR(IR)*IPT`
3. 占据数检查 `ICHKQ1(IC, IR)`

对应代码位置：

- [hfsgg_mpi.f90](../../src/appl/rhfs90/hfsgg_mpi.f90)

这样做的目的，是避免大量本来会在 `ONEPARTICLEJJ` 内部立刻返回空结果的调用。

### 4.2 保留轻量诊断计数，删除高频重计时

第一轮 profiling 里曾加入更细粒度的 `SYSTEM_CLOCK` 统计，虽然帮助定位了瓶颈，但也引入了明显的测量开销。

因此在正式优化版本中：

- 保留对筛选路径有诊断价值的计数器
- 去掉高频细粒度计时

该阶段保留的关键计数包括（后来已在 `3d067e7` 删除）：

- `cand`
- `active`
- `ptrig`
- `ppar`
- `pocc`
- `ia0`
- `diag`
- `offd`
- `zia0`
- `zelem`

这些计数器用于确认工作项是在何处被筛掉，而不会像高频计时那样明显拉高运行时间。

## 5. 优化效果

### 5.1 总时间变化

同一类输入的代表性结果如下。

优化前：

- `HFSGG_MPI kernel wall time (max rank): 2048.128 s`
- `HFSGG_MPI reduce wall time (max rank): 15.238 s`

优化后：

- `HFSGG_MPI kernel wall time (max rank): 1429.850 s`
- `HFSGG_MPI reduce wall time (max rank): 35.083 s`

结论：

- `kernel` 从 `2048.128 s` 降到 `1429.850 s`
- 绝对减少 `618.278 s`
- 降幅约 `30.2%`

虽然 `reduce` 时间有所波动，但总瓶颈仍然是本地内核计算，这次优化的主要收益来自计算侧而不是通信侧。

### 5.2 `ONEPARTICLEJJ` 调用量变化

以单个 rank 的典型数据为例，优化前后对比：

- 优化前：`active ~= 1.51e9`
- 优化后：`active ~= 4.00e8`

也就是说，真正进入 `ONEPARTICLEJJ` 的调用数下降到了原来的约四分之一。

对应地，每个 rank 的 `onep` 时间也从大约：

- `1570 - 1590 s`

下降到：

- `385 - 425 s`

这说明前置筛选直接压缩了最昂贵的 `ONEPARTICLEJJ` 调用开销。

### 5.3 筛选路径变化

优化后日志显示：

- `ptrig = 0`
- `ppar = 0`
- `pocc ~= 1.11e9`
- `ia0 = 0`

这说明：

- 三角条件和宇称在当前这组问题上不是主要过滤来源
- 占据数条件 `ICHKQ1` 是最关键的前置筛选
- 原本大量 `IA = 0` 的空返回已被外层 `ICHKQ1` 成功提前截断

## 6. 当前剩余瓶颈

本次优化之后，新的主要浪费已经不再是 `IA = 0` 的空调用，而是：

- `active ~= 4.0e8`
- `zelem ~= 3.7e8`
- `nzacc ~= 2.8e7`

这说明当前进入 `ONEPARTICLEJJ` 的项中，仍然有大约九成以上在 `TSHELL -> ELEMNT / ELEMNTGJ / ELEMNTDGJ` 这一阶段被判零，最后没有进入最终累加。

因此当前的主瓶颈已经从：

- “空的 `ONEPARTICLEJJ` 调用”

转移为：

- “`ONEPARTICLEJJ` 返回非空壳层信息，但后续元素构造仍为零”

## 7. 后续优化建议

本节是第一阶段建议；后续 `a2ec419` 已加入预计算和零元素筛选。
当前保留项见[分支总结](1.0.1-dev.1_optimization_summary.md)及当前源码，
不要把下面建议全部当成未实施任务。

如果继续做下一轮优化，优先级建议如下：

1. 研究能否在 `ONEPARTICLEJJ` 返回后、进入 `TSHELL -> ELEMNT` 累加前，再加入一个更便宜的零贡献判断。
2. 如果无法进一步前移筛选，则应直接优化 `TSHELL` 到 `ELEMNT/ELEMNTGJ/ELEMNTDGJ` 的构造过程。
3. 在后续 profiling 中，尽量避免在最深热路径中加入过于频繁的 `SYSTEM_CLOCK` 调用，以免测量本身干扰性能判断。

## 8. 小结

本次优化的核心不是改 MPI 通信，而是把原本位于 `ONEPARTICLEJJ` 内部的廉价占据数筛选前移到 `HFSGG_MPI` 外层，减少了大量无效的单电子矩阵元调用。

从实际结果看：

- `HFSGG_MPI kernel` 降低约 `30.2%`
- `ONEPARTICLEJJ` 相关时间降低约 `75%`
- `IA = 0` 类型的空调用基本被消除

因此，这次测量表明：在该样例的 `RHFS_MPI` 的 `HFSGG_MPI` 路径中，减少无效工作项带来了主要提速；这一结果不代替当前版本的重新测量。
