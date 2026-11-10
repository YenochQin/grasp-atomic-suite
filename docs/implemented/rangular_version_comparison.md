# 原版与新版 rangular 对照运行脚本

日期：2026-11-10。入口：[compare_rangular_versions.py](../../scripts/compare_rangular_versions.py)。
范围：**普通 rangular 的串行原版与串行/OpenMP 候选**。脚本接受两个已构建可执行文件，不自动编译，不运行 MPI 或提交 SLURM。
新版是已接入 `--csf-descriptors` / `--timings` 的 suite 版本。

## 1. 最小用法：原版文本与新版文本

从 suite 根目录执行；把示例中的可执行文件和数据路径替换为实际路径：

```fish
../graspkit-tools/.venv/bin/python scripts/compare_rangular_versions.py \
  --reference /path/to/original/rangular \
  --candidate /path/to/new/rangular \
  --csf /path/to/rcsf.inp \
  --output build-rangular-comparison/text-case \
  --repeat 3
```

默认给两版同样的 `y` 输入，表示完整相互作用。原版不接收任何新增命令行选项；
两版分别从自己的工作目录读取快照 `rcsf.inp`，新版另外启用 `--timings timings.csv`。
脚本不修改原数据或可执行文件。结果目录必须尚不存在；再次运行需选择新目录。

依赖为 Python >=3.11、NumPy 和同目录的 `mcp_compare.py`。
仅比较文本或使用预制二进制时，运行器本身不导入 PyArrow。

## 2. 同时测试新版描述符输入

已有 rCSFs V2 Parquet 和匹配 header 时：

```fish
../graspkit-tools/.venv/bin/python scripts/compare_rangular_versions.py \
  --reference /path/to/original/rangular \
  --candidate /path/to/new/rangular \
  --csf /path/to/original.c \
  --descriptors /path/to/descriptors.parquet \
  --header /path/to/source_header.toml \
  --verifier /path/to/new/csf_descriptor_verify \
  --output build-rangular-comparison/descriptor-case \
  --repeat 3
```

这会依次准备并核对输入，然后在每轮分别运行：原版文本、新版文本、新版二进制。
Parquet 预处理只执行一次，额外需要 PyArrow；其耗时单列，不计入 rangular 墙钟。
若已有交换文件，用 `--binary /path/to/csfs.bin` 替代 `--descriptors ... --header ...`。
`--verifier` 默认在新版 rangular 同一目录寻找 `csf_descriptor_verify`。

描述符模式先逐项核对其内部数组与原文本；验证器的返回码和最终成功标记都必须通过。
二进制运行时不放置 `rcsf.inp`，结束后才将已验证的文本快照加入该目录，供 MCP 对比建立基组身份。
这个动作在 `execution.json` 中标为 `text_provenance_added_after_execution`。
交换格式及支持态范围见[描述符输入说明](rangular_descriptor_input.md)。

使用 [rangular_omp](rangular_openmp.md) 时，将 `--candidate` 指向它并增加 `--candidate-threads 8`。
该选项同时设置候选的 OpenMP 环境和 `--threads` 参数；原版及数学库仍单线程。
省略此选项时，运行器保持原有单线程默认，不能仅设置父 shell 的 OMP_NUM_THREADS 就认定候选已多线程运行。

## 3. 一阶相互作用、超时与判据

将完整交互输入放入文本文件。例如三个 block 各使用 `ICCUT=1`：

```text
n
1
1
1
```

通过 `--stdin-file /path/to/rangular.stdin` 指定。双方必须使用兼容的交互协议；
如果更老的原版会询问其他设置，需要提供适合双方的输入，不能据相同输入文件推断解释必然相同。
完整 transcript、输出头部及实际比较结果都会保留。

- 默认 `--atol 0 --rtol 0`，按 CSF 对、积分标签和每项贡献比较，保留重复贡献数量。
- 明确需要容差时可指定 `--atol` / `--rtol`；应在实验前约定，脚本不会根据结果自动放宽。
- 缺失贡献、稀疏结构或元数据差异仍是失败；不能用数值容差忽略。
- `--keyorb` 默认 215，必须与双方构建的积分标签基数一致。文件本身不保存该常量。
- `--timeout SECONDS` 限制每次子进程（包括数组验证）的时间；默认不限时。
  超时/中断保存日志，POSIX 上终止本次子进程组，不留下后台计算。

退出码：`0` 为所有对照通过，`1` 为有效输出存在差异，`2` 为运行/输入/格式错误，`130` 为用户中断。
旧 Fortran 提前 `STOP` 即使返回零，只要 MCP 不完整也会失败。

## 4. 报告与计时

```text
结果目录/
  summary.json                         总结、身份、环境、各次结果和墙钟统计
  inputs/                              原 CSF、交互输入、可选描述符/header/二进制快照
  array_verification/                  可选的内部数组验证日志
  trial_001_reference_text/            原版 stdout/stderr、MCP、execution.json
  trial_001_candidate_text/            新版文本输出，加 timings.csv
  trial_001_candidate_descriptors/     可选新版二进制输出，加 timings.csv
  comparison_001_candidate_text.json   与当轮原版的逐项比较
  comparison_001_candidate_descriptors.json
  reference_repeatability_002.json     后续轮次还检查原版输出是否可重复
```

报告保存两个可执行文件及输入 SHA-256，并在每次计算前后检查可执行文件是否变化。
若双方实际为同一文件内容，`same_executable_bytes=true` 明确记录，不把它当成独立旧/新版本证据。
比较同时报告语义判据与 `byte_identical`；记录布局不同但数值相等可通过语义判据，字节一致性单列。

外部墙钟仅涵盖子进程运行；输入复制、预处理、数组验证、MCP 校验与比较均在该区间之外。
重复运行按轮次轮换执行顺序，汇总中位数。只有所有对照通过才报告
`reference_over_candidate_wall_ratio`，含义为原版墙钟中位数除以对应新版中位数。
新版阶段时间保存在各次 CSV 及 `stage_wall_seconds`，整体区间包含子阶段，不能再次求和。

默认将 OpenMP/OpenBLAS/MKL/BLIS/Accelerate/NumExpr 的线程环境设为 1；
`--candidate-threads` 仅覆盖候选的 OpenMP 线程数，实际每次环境记录在 execution JSON 中。
不修改绑定、CPU 频率或系统缓存。输入复制和验证可能预热页缓存，因此这些结果不是冷缓存测试。
建议使用相同工具链、优化等级、容量及科学设置；尤其不要把旧 Debug 与新 Release 的比值当成改造收益。
大量 MCP 文件的校验可能耗时，但不会计入可执行程序墙钟。

脚本执行的是角系数与文件接口回归，不证明完整 SCF/RCI、活性空间收敛或物理精度。

## 5. 本轮验证

[七项维护测试](../../test/test_compare_rangular_versions.py)覆盖真实三路径/重复运行、预制二进制与一阶输入、
系数差异、零退出但缺输出、超时、验证器旧 STOP、非法参数和已有结果目录保护。
临时路径含空格，命令采用参数数组传递；不通过 shell 拼接文件路径。

另外使用本机已有旧 GNU 15 Debug 二进制与本次新版 Check 二进制，对 `test/jj2lsj/f3.c`
及其真实 rCSFs 描述符执行两轮三路径对照：四组旧/新比较均数值和 MCP 字节一致，原版跨轮次也一致。
此旧二进制是现有构建产物，身份以结果中的 SHA-256 为准，本轮没有重新建立其完整源码 provenance。
记录在忽略的 `build-descriptor-check/audit/version-comparison-f3/`。

同一旧 Debug 对新版 Release 的另外一次运行，默认零容差正确报告了每组 274 项末位差异，
最大绝对差约 `4.996e-16`，没有缺失贡献；新版文本和描述符路径表现相同。
未通过放宽容差将该结果改记为通过，也不据这些小样本墙钟宣称性能提升。

复现维护测试：

```fish
../graspkit-tools/.venv/bin/python test/test_compare_rangular_versions.py \
  --rangular build-descriptor-check/bin/rangular \
  --verifier build-descriptor-check/bin/csf_descriptor_verify -v
ctest --test-dir build-descriptor-check -R '^rangular_version_comparison$' --output-on-failure
```
