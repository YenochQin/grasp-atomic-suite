# jj2lsj 分阶段计时与 ASF 排序优化

状态：2026-11-04 已实现分阶段计时、稳定归并排序和 ASF 级排序并行。
修改前基线：`9807b2806e9821e5d383b5de446e68f3d697d5c9`。
跨 block 异步流水线经过真实输入计时评估后暂缓；唯一标签和输出保持有序串行。
这属于程序兼容性和性能工作，不是原子物理收敛验证。

## 使用方式

仍使用 `jj2lsj_omp`，问答、文件格式和 Tools 生成器接口不变。
默认在转换结束后对各 ASF 的 LS 系数进行稳定排序，多个 ASF 可并行排序。
只有一个 ASF 或没有 OpenMP 时自动顺序执行。`OMP_NUM_THREADS` 控制线程预算。

```sh
# 在已有输入目录运行；替换安装路径和问答文件。
OMP_NUM_THREADS=8 OMP_DYNAMIC=FALSE JJ2LSJ_PROFILE=1 \
    /path/to/grasp-atomic-suite/bin/jj2lsj_omp < jj2lsj.inp > jj2lsj.log
```

`JJ2LSJ_PROFILE=1` 在当前目录创建或覆盖 `jj2lsj_omp.profile.csv`。
默认不记录计时；关闭开关时不会更新或删除旧计时文件。
同目录连续转换多个计算时，应在下一次调用前保存该文件。

`JJ2LSJ_POST_SERIAL=1` 仅让新排序按 ASF 顺序执行，可用于隔离排序并发的影响；
`JJ2LSJ_SERIAL=1` 同时关闭转换和排序的并发。两种模式都使用新的稳定归并排序，
不是恢复旧的二次复杂度插入排序。与旧算法对照应使用改动前的可执行文件。

## 计时内容及解释

计时模块为 [jj2lsj_profile.f90](../../src/appl/jj2lsj90/jj2lsj_profile.f90)，使用 `SYSTEM_CLOCK`。
CSV 列为 `phase,block,level,items,seconds`；`block`、`level` 为 0 表示不适用。
能级编号为输入 mixing 数据中的全局编号。`total` 包括交互等待，应使用固定问答文件做性能对照。

| phase | 范围或 items 含义 |
| --- | --- |
| `csf_read_parse` | CSF 文件读取和解析；items 为 jj-CSF 数 |
| `jcoup` | 耦合表分配及建立；items 为 jj-CSF 数 |
| `mix_read` | 混合系数加载；items 为输入能级数 |
| `screen`、`ls_basis` | 筛选和 LS 基生成；items 分别为入选 jj-CSF 数、LS-CSF 数 |
| `workspace`、`release` | 主要转换工作区分配、释放 |
| `jj_rank` | 单个 ASF 的 jj 主要成分选择 |
| `transform` | 一个 block 的全部 ASF 的 jj–LS 转换 |
| `ls_rank_wall` | 整个 block 的新排序墙钟时间，包括线程启动与等待；items 为请求的线程组大小 |
| `ls_rank_asf` | 各 ASF 的排序及原顺序范数累加时间；items 为 LS-CSF 数 |
| `ls_rank_gather` | 把排序索引整理为已有输出流程需要的数组 |
| `candidate_select` | 打印候选选择，包含范数警告检查；items 为打印成分数 |
| `jj_stdout`、`ls_stdout` | 控制台格式化与输出，含主要成分描述准备 |
| `label_format`、`label_write` | `.lsj.lbl` 的标签字符串准备，以及其余记录格式化/写入 |
| `csf_write`、`vector_write` | 可选 `.lsj.c` 和 `.lsj.j` 输出 |
| `unique_read` | 唯一标签各组的文件读回；items 为组内 ASF 数 |
| `unique_assign_summary` | 唯一标签协调及 `.uni.lsj.sum` 写入；未将每条汇总写入单独计时 |
| `unique_write` | `.uni.lsj.lbl` 的最终输出 |
| `unique_total` | 整个唯一标签过程，包含上述三个子阶段及分配等开销 |
| `output_flush` | 计时模式中测量输出大小前的缓冲刷新 |
| `*_bytes` | 已生成的对应输出文件大小，items 为字节数、seconds 为 0 |
| `total` | 主程序完整墙钟时间，包含诊断开销 |

这些时间不是一组可全部相加的互斥区间：`ls_rank_asf` 在并行情况下重叠，
`unique_total` 包含其子阶段，`total` 包含全部工作。
格式化与写入时间包含运行库开销，不能直接当作磁盘设备服务时间。
原排序计时基线的 `ls_sort` 对应优化后的 `ls_rank_wall` 加 `ls_rank_gather`；
不能把所有 ASF 的并发时间相加后当作新的墙钟时间。

## 数值和线程契约

[jj2lsj_postprocess.f90](../../src/appl/jj2lsj90/jj2lsj_postprocess.f90) 对 LS 索引进行迭代归并排序：

- 按系数平方降序排列，只有右项严格大于左项时才先取右项，保留相等权重的原始 LS 顺序。
- 正负系数和正负零均从原向量按索引读取，不对原向量进行重排或筛选。
- 范数按原始 LS 顺序逐项累加，保持原串行求和次序。
- 每个 ASF 有独立的排序索引和临时数组，复杂度为 `O(N_LS log N_LS)`。
- 非有限输入系数会明确报错；不对 NaN/Inf 的未定义排序行为作兼容承诺。

原先的 `ithresh`、LS 枚举、转换核、共同提前停止点和打印阈值保持原算法。
转换结束后才开始 ASF 排序；排序完成后按原能级顺序执行格式化和输出，
最后仍由原 `uniquelsj` 从 `.lsj.lbl` 读回已舍入的数值来协调唯一标签。
没有绕过原文件精度，也没有让各 ASF 独立抢占唯一标签。

并行排序和转换不重叠，也不建立嵌套线程组。多 ASF 排序沿用转换的最大线程组大小，
避免主动改变相邻活动并行区的团队大小；有几个 ASF 才有几个排序任务，闲置线程不会增加有效并行度。
继续使用 `OMP_DYNAMIC=FALSE`。线程私有数据跨并行区的保留条件见
[OpenMP threadprivate 规范](https://openmp.org/spec-html/5.2/openmpse24.html)。

额外持久排序存储约为 `sizeof(integer) × N_LS × 所选ASF数`，
每个活动排序任务还需一份实数平方数组和整数临时索引。工作区在当前 block 结束时释放。
未改变原 `EVEC` 的全局布局，也没有同时保留多个 block 的结果。

## 验证及范围决定

验证包含稳定排序对旧插入排序的独立对照、相同权重、正负零、原顺序范数的位模式、
84,000 行零系数平台、非有限值错误路径；现有系数/文件回归进一步检查多 ASF 并发、
计时开关、输出字节数、唯一标签、`.lsj.T` 读写和早停。
另用改动前程序对照 15 组小型输入/模式的结果文件，包括完整向量和变换文件。
Release、带 `-fcheck=all` 的 Debug、无 OpenMP 构建均进行相应检查。

真实输入的设置、结果与新增方案的取舍见
[2026-11-04 后处理验证报告](../reviews/jj2lsj_postprocessing_validation_2026-11-04.md)。
基线实测中唯一标签和标签输出只占极少时间，稳定排序之后可供后台隐藏的工作更少；
因此本次没有加入跨 block 异步流水线，也没有并发调用仍带静态工作区的格式化例程。
如果后续多 ASF 输入或不同文件系统显示这些阶段确有显著开销，再按
[原方案](../plans/jj2lsj_postprocessing_parallelization_2026-11-03.md)评估独立上下文、有界缓冲和异步提交。
