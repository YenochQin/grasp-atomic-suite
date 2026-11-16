# 普通 GRASP 构建范围与编译选项

更新：2026-11-12，基于 `46d1d54` 的工作树增量。按当前工作范围，CSFG 隔离为可选扩展；生产选项和诊断选项分开。
本页取代旧迁移计划中保留 legacy/隐式静态存储的构建建议。

## 范围

`GRASP_BUILD_CSFG` 默认 **OFF**。CSFG 专用源代码保留在已有独立目录，唯一构建入口集中到
[CSFGExtension.cmake](../../cmake/CSFGExtension.cmake)，默认不编译其程序、库或源码测试。
普通构建中的检查会拒绝意外进入编译命令清单的 CSFG 源文件。
本轮不继续实现或验收 CSFG；旧 CSFG 结果属于历史证据，不代表当前隔离分支的新验证。

## 三种配置

| 配置 | 构建类型 | 用途与主要选项 |
| --- | --- | --- |
| `f2018-release`，默认 | Release | 生产计算：`-O2`；并行目标另加 `-fopenmp`。不附带语言诊断、运行时检查或警告采集 launcher。 |
| `f2018-check` | Debug | 标准/接口/边界等验证：`-std=f2018 -pedantic-errors -Wall -Wextra -Werror=implicit-interface -Werror=implicit-procedure -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow`。 |
| `f2018-audit` | Release | 优化后诊断：上述标准/接口/警告选项配 `-O2`；无隐式 SAVE、无运行时检查。 |

源码和所有配置仍要求 GNU Fortran 15。`legacy` 配置已移除，所有配置拒绝
`-std=legacy`、`-fno-automatic`、`-fallow-argument-mismatch`、`-fall-intrinsics` 及全局默认 kind 扩宽。
没有引入 `-ffast-math` 或 `-Ofast` 等会改变数值假设的选项。

已删除冗余的 `-fno-fast-math`、`-ffp-contract=off`、为旧存储配置兜底的 `-fautomatic`，
以及 OpenMP 已隐含的 `-frecursive`。需要持久的状态在源码中明确表达；例如 SIXJ 的 RACA 缓存
显式 SAVE，在并发对象中仍为 threadprivate，没有取消缓存。

`-I/-J/-c/-o` 等编译/模块路径参数仍是正常构建机制。
数值测试目标可以拥有专门的检查/poison 选项，它们不传给生产计算目标。
CMake 的范围核查与测试设施需要 Python；生成的 Fortran 数值程序不链接 Python 或 Rust。

## 使用

`./configure.sh` 默认配置全部保留的非 CSFG 程序：普通串行目标（包括
`rwfnestimate`）、MPI 目标及两个 OpenMP 目标。随后执行
`cmake --build build --parallel 4`，可执行文件统一输出到 `build/bin`。
同一源码与编译器的已有 `build` 会更新配置，因此旧缓存中的串行关闭状态也会被纠正。
若需要精简范围，显式传入 `--no-serial`、`--no-mpi` 或 `--no-openmp`，每次重配均需保留
这些选项。`--serial` 继续兼容；直接调用 CMake 时仍需设置
`-DGRASP_BUILD_SERIAL_APPS=ON`。

生产 OpenMP 构建，使用新目录：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC $FC
./configure.sh --serial --no-mpi --no-csfg \
    --profile f2018-release --build-dir build-rangular-production \
    -DGRASP_JJ2LSJ_OPENMP=OFF -DGRASP_RANGULAR_OPENMP=ON
cmake --build build-rangular-production \
    --target rangular rangular_omp csf_descriptor_verify -j4
```

检查配置可用 `--debug --profile f2018-check`，优化诊断配置用 `--profile f2018-audit`。
旧目录中缓存的 profile/警告采集设置不会被当作新配置的默认值；尤其旧 Release 曾启用的
`GRASP_CAPTURE_WARNINGS=ON` 不属于当前生产配置。优先用独立目录。

普通验证矩阵独立覆盖 serial、MPI、jj2lsj OpenMP 和 rangular OpenMP。
serial 开启时，每种 profile 有 MPI 与两个 OpenMP 开关的 8 种组合；serial 关闭时
rangular 目标不存在，只保留 MPI 与 jj2lsj OpenMP 的 4 种组合。三种 profile 共 36 种配置；
CSFG 始终关闭。24 个 Check/Audit 组合均须有各自的警告预算，回归测试检查基线覆盖。
生产结果明确标记未采集警告，不能以生产构建
没有诊断输出宣称告警已清零。矩阵工具仍需实际运行全部组合后才能报告完整矩阵验收。

独立开关补齐保留了已有 18 个预算，新增 6 个配置由本机 GNU 15.3.0 实际编译采集，
记录在基线的 `independent_openmp_budget_update`。GNU 15.3 在 macOS arm64 上将
`rci/qed.f90` 和 `rci_mpi/qed.f90` 中既有的 K 未初始化警告定位到使用处，
Linux 基线定位在声明处；新增 Audit 预算为同一警告保留这两种确切身份，各限一次。
源文件哈希未变，其他诊断没有超过已有身份的最大预算。新增配置同时保留原 Linux
对应目标的预算，包括本机未发出的 BIOTR1 告警，避免用单平台采集结果清除既有风险；
这些 QED 和 BIOTR1 风险仍待单独复核。

## 本次 flag 对照

当前环境的 GNU 15 默认目标为通用 x86-64，FMA 默认关闭。对同一 23,459-CSF 二进制输入，
用相同的 8 个线程交替运行整理前的构建与精简生产构建，每种三次。
目标是确认 MCP 完全一致并观察耗时，不预设精简会加速。

| 配置 | 三次墙钟 | 中位数 |
| --- | --- | ---: |
| 整理前 | 24.945 / 25.910 / 25.274 s | 25.274 s |
| 精简生产 flags | 24.437 / 25.464 / 25.662 s | 25.464 s |

六次 MCP 与同一基线逐字节一致。两组范围重叠，中位数差约 0.75%，未体现保留额外 flags 的稳定优势；
不把这个小差异当作可靠加速或退化结论。按当前取向采用精简生产配置。
记录在忽略的 `build-rangular-production/audit/flags-23459-v1/`。
不同 ISA/编译器默认值仍需重新验证数值与性能；没有据此启用本机专属指令或激进浮点选项。

当前选定配置的回归：普通 Check/Audit/生产 Release 各 44 项，关闭 OpenMP 的 Check 42 项，
普通 MPI 生产配置 78 项，均通过。这些是独立开关补齐前的记录，未运行完整 36 组合矩阵。

## 独立开关的本机验证（2026-11-11）

使用 macOS arm64、GNU Fortran 15.3.0、OpenBLAS 和共享 Python 环境，实际编译补齐的
6 个配置，最终基线的警告检查全部通过。CTest 结果如下，两个诊断 profile 分别运行：

| jj2lsj OpenMP | rangular OpenMP | MPI | Check | Audit |
| --- | --- | --- | --- | --- |
| ON | OFF | OFF | 42 通过，无跳过 | 42 通过，无跳过 |
| ON | OFF | ON | 72 通过，4 跳过 | 72 通过，4 跳过 |
| OFF | ON | ON | 74 通过，4 跳过 | 74 通过，4 跳过 |

合计 376 次通过、16 次跳过，没有执行失败。每个 MPI 配置的跳过项相同：
`p0_trace_mpi_1/2`、`strict_mem_restart_mpi_1/2`，其上游 memory solver 需要 Linux
`/proc/meminfo`。严格验收仍将跳过判为未通过；没有更改测试或放宽此规则。
上述四项属于随后已否决并移除的 strict/diag 变体，现不再注册；该段仅保留当时的验证结果，
不能用它证明当前配置验收。现行范围见[移除说明](rejected_rmcdhf_variants.md)。
可追溯汇总为忽略目录 `build-review-openmp-independent/local-validation-summary.json`。
