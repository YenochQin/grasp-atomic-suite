# rangular 保守角向系数复用

日期：2026-11-14。基于 `3f82922d2caa4de061c7f62c1b604770d0092847` 的工作树实现。
范围：普通串行 `rangular`，文本与 rCSFs 描述符入口共用；默认关闭。
MPI 不启用此缓存；rangular OpenMP 与 CSFG 已移除。实现入口为
[rangular_reuse.f90](../../src/apps/rangular/rangular_reuse.f90)，
接入点为 [mcp_gg.f90](../../src/apps/rangular/compute_mcp_angular_coefficients.f90)。

## 使用

在新的计算目录内运行已构建的串行程序，交互输入与原版相同：

```fish
/path/to/rangular --csf-descriptors /path/to/csfs.bin \
  --reuse-cache-mib 64 --reuse-stats reuse.csv --timings timings.csv
```

文本输入可用 `--csf-text /path/to/rcsf.inp`，或使用工作目录的 `rcsf.inp`。
`--reuse-cache-mib` 接受 0–4096，0 保留原计算路径；正值启用复用。
统计与计时文件必须尚不存在。MPI 保留原入口，不接受这些选项。

## 精确等价条件

缓存按 CSF 对建立，不能只按某个 CSF 的 κ 序列建立。两个输入都已转换成相同的
`IQA`、`JQSA`、`JCUPA` 后才进入缓存，不改变 rCSFs 描述符格式或补全规则。

1. 按轨道表原顺序构造 bra/ket 的非空轨道并集，仅去掉两侧同时为空的轨道。
   每个位置保存 κ、双方占据数，以及双方完整子壳层状态 `(v,w,2J+1)`。
   同一并集位置在双方是否占据，明确保留共同轨道身份。
2. 包含对角/非对角标志，并严格比较双方完整 `JCUPA`，包括最终 J/宇称。
   不交换 bra/ket，不按 κ 排序，不重新耦合，不尝试推断相位。
3. 主量子数 n 不进入键。查询使用哈希加开放寻址，但哈希相等并不足以复用：
   必须逐项比较完整轨道键和原始耦合数组。耦合哈希也只用于索引。
4. 未命中调用原 `ONESCALAR` / `RKCO_GG(...,CORD,0,1)`，保存未经 MCP 后处理的
   T、V 系数序列及抽象轨道编号；命中只恢复浮点值并映射回当前真实轨道编号。
   保持原贡献顺序和重复项，不对系数求和或重新计算浮点运算。
5. 调试输出启用、分配失败、预算不足，或内核返回无法映射的轨道时回退原计算。
   不补全输入中缺失的 CSF，不把同族 CSF 合并为一个混合系数。

现有 MCP 截断、F0 排除、轨道标签规范化、稀疏结构、排序与写出流程保持原样。
因此输出仍是普通展开的 `mcp.*`，大小不会因为缓存而变小；下游求解器无需修改。
所有原有 CSF 对仍被访问，复杂度没有变成 CSFG 的生成元对遍历。

## 内存与统计

每个对称性块重新建立缓存，处理完立即释放。缓存用固定容量的桶、条目、键和
系数数组；没有逐对动态分配。填满任一数组后停止插入，已存条目仍可命中，未命中
继续调用原内核；不驱逐旧条目，也不把全部 CSF 对索引放入内存。

`reuse.csv` 每块记录：

| 字段 | 含义 |
| --- | --- |
| `pairs` | 实际访问的 CSF 对数，已应用原 DIAG/LFORDR 筛选 |
| `hits` / `misses` | 查询成功 / 需要原内核计算的次数；无系数结果也可以缓存 |
| `bypasses` | 缓存不可用、直接使用原内核的次数 |
| `not_stored` | 未命中后因容量或映射限制未保存的次数 |
| `entries` | 已保存的精确键数量，不是整个块的不同键总数 |
| `probes` | 查询检查过的非空哈希桶数量 |
| `cache_reserved_bytes` | 本块为缓存数组实际分配的存储字节数 |
| `cache_used_bytes` | 桶和耦合哈希，加已使用条目/键/系数的字节数 |
| `cache_limit_bytes` | 用户指定预算 |

复用率定义为 `sum(hits)/sum(pairs)`，不能把条目数直接当成整个计算的唯一模式数。
预算包含每个 CSF 的耦合哈希和桶，但不包含原有 CSF/MCP 数组、Fortran 分配器开销、
固定大小局部工作数组或运行库。它是缓存数组预算，不是整个进程的内存上限。

[现有对照程序](../../scripts/compare_rangular_versions.py) 新增
`--candidate-reuse-cache-mib`；每次候选运行保存并校验 `reuse.csv`，汇总到 JSON。
`--reference-timings` 允许支持计时参数的 suite 基准同时记录阶段时间。
Linux/macOS 使用针对当前子进程的 `wait4.ru_maxrss` 测量进程峰值 RSS，统一为字节，
不使用累计的子进程最大值；其他平台报告 unavailable。总墙钟包含最多约 10 ms 的
退出轮询延迟，小于此量级的样本不能用于判断加速。

## 对照与适用边界

```fish
../graspkit-tools/.venv/bin/python scripts/compare_rangular_versions.py \
  --reference /path/to/baseline/rangular --candidate /path/to/candidate/rangular \
  --csf /path/to/rcsf.inp --binary /path/to/csfs.bin \
  --verifier /path/to/csf_descriptor_verify \
  --candidate-reuse-cache-mib 64 --reference-timings \
  --repeat 3 --output build-angular-reuse/comparison
```

默认零容差，另外检查文件是否字节一致。为隔离复用开关，也可把 reference 和
candidate 指向同一个二进制；JSON 会明确记录这一点。基准若是早期上游程序，
省略它不支持的 `--reference-timings`。

[专项回归](../../test/rangular/test_reuse.py) 覆盖同 κ 不同轨道、共享轨道/无共享轨道、
不同耦合链与子壳层项、交错轨道顺序、不完整轨道族、相同 J/宇称的相邻块、描述符/
文本输入、一阶相互作用、缓存耗尽、非法参数和已有文件保护。CTest 名为
`rangular_angular_reuse`，GNU 15 Release 与 Check 均验证。

收益取决于允许保持顺序的精确等价模式数量、遍历顺序、缓存容量及原内核耗时。
没有径向同族结构时，哈希查询会增加开销，填满缓存也可能错过后来出现的可复用族。
此版本保留显式开关，不能据“可无损复用”推断任意 CSF 数据都会加速。
性能实测与证据见[本轮报告](../reviews/rangular_reuse_benchmark_2026-11-14.md)。
