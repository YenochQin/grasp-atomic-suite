# 原生描述符迁移验收记录

2026-11-17 的原生输入迁移以 Suite `5806884` 和 rCSFs `fea061f` 为源码基线，验证对象为基线之后的迁移增量。
rCSFs 实现随后记录为提交 `f17a862d77c45164a6d31ae1e9cfeb4efbd04cf6`。
2026-11-28 的 V3 复验见本文末节。[现行输入契约与构建方式](../implemented/design-0001-native-csf-descriptors.md)

## 范围和对照

本任务验证输入表示、数组、应用状态、MPI 协议和小型数值流程的一致性。
未做生产规模性能、活性空间收敛或最终物理精度研究。

保留程序均已接入：串行/MPI rangular、四个 RMCDHF、串行/MPI RCI、biotransform、
transition/phase、JJ2LSJ 串行/OpenMP、gj/rhfs/ris4/rdensity 串行/MPI，以及 header-only rwfnestimate。
旧文本 reader 仅在不安装的 `test/descriptor_oracle` 中使用；生产库无 CSF 文本读取入口。

本地忽略目录 `build-rust-descriptors/` 保存构建、CTest、源码身份及对照产物。
基线来自 `git archive HEAD` 的独立源码副本，使用同一 GNU Fortran 15 工具链构建，不覆盖安装目录。
测试输入、原生 producer、程序和产物的身份记录属于本次回归证据，不等于给旧 MCP 补充了协议级身份绑定。

## 有针对性的验证

- Rust：必需 footer 键的缺失/重复/错误值，真实 Int8 schema 与空值检查，header 字节绑定、块计数、
  截断、显式闭核、批次/row group 顺序、相邻同对称性块、独立句柄，以及 C ABI 的缓冲容量和关闭/重开。
- Fortran：保留完整三数组与历史合法文本 oracle 的逐项对照；仅占据入口与完整入口精确一致，
  深层物理错误返回、关闭和重开；纯闭核及没有 Python/Rust 命令可用的运行路径。
- RMCDHF：用实际 rangular 产出 MCP，再调用各版本生产 SETMCP/CSF 加载入口。
  四版本占据字节精确匹配，块数组下界保持 0、前缀/长度未覆盖，JPgg 符号与值正确，
  无全基组 JQSA/JCUPA 分配。MPI 1/2/4 ranks 覆盖 root 输入错误和 worker MCP 元数据冲突。
- EAL 标准权重：两个 J 块、三个 CSF 的解析期望权重为 `[0.2,0.2,0.6]`，
  四个版本均在耦合指针为空时得到该结果。
- MPI 公共层：1/2/4 ranks 强制每次只传一列；worker 负错误码、不同形状明确非零失败，设置超时防止挂起。
- 双基组：左右轨道表不同、多个轨道插入、块前缀、量子数映射和次序冲突；不生成 `SLASK`。
  保持两个输入的相对次序，不能在插入多个轨道后改变某一态的耦合顺序。
- 端到端：小型一电子 RCI、biotransform、transition/phase 与文本基线比较；
  `.cm`、`.bm`、`.bw` 和跃迁输出按完整产物核对，不只看退出码。
- 固定一次迭代的 He 软件回归：两个串行版本、两个 MPI 版本的 1/2/4-rank 运行，与相同程序/rank 的
  基线径向输出逐字节一致，均正常到达设定的一次迭代上限。该结果不是 SCF 收敛结论。
  追溯目录为 `scf-audit-v2/`、`scf-mpi-audit/`；后者的 `summary.json` 记录各版本/rank 的输出哈希。

RCI 原数值内核要求每个块的 CSF 数不少于 rank 数；保留这个限制，并使不支持的配置非零失败。
多 rank RCI 端到端夹具使用满足该约束的块，不通过删除内核保护来制造通过结果。

## 回归发现的既有问题

以下问题在输入迁移前的源码/独立基线中同样存在。为完成小型 Check 回归，
对照副本与候选实现采用相同修正；不能把这些修正的效果归因于 Parquet 输入：

| 问题 | 修正 |
| --- | --- |
| biotransform 的 `csf_mixing_coefficients(1:500)` 对小基组越界 | 传入实际数组范围，CITRAG 的维度仍由 NCSF/NCIV 决定 |
| RCI `QED` 把按轨道填写的 UCF 声明为 `ucf(1)` | 实现和显式接口同步改为 assumed-shape，计算表达式不变 |
| MPI biotransform 的 getmix 要求 128 字符，但调用实参更短 | 实现和接口采用实参长度 |
| MPI biotransform 空 rank 在 CITRAG 前跳过循环 | 空 rank 继续参与已有的零贡献 collective，避免通信不匹配 |
| MPI transition 的首次 strsum 使用未初始化 ilbl | 明确选择正常输出模式 0，与串行入口一致 |

直接数组合并还采用共享的稳定合并算法，避免旧插入循环在多次插入后漏过已移动的公共轨道。
单独的次序/映射测试覆盖该路径，未引入轨道或 CSF 的任意排序。

## 告警门禁边界

未增加旧告警预算、关闭诊断选项或放宽编译器检查。
现有 `gfortran15-warning-baseline.json` 在迁移前的独立 HEAD 构建上已报 source coverage changed：
其该配置预算为 1,897 条、1,873 个文件，而当前基线源码的实际编译覆盖和诊断数量已经不同。
候选的旧 reader 退役及新增测试也改变覆盖范围。因此这里不宣称该全局历史告警门禁通过。

另修正了检查脚本漏列现有 `test.jj2lsj_pointer_defaults` OpenMP 测试目标的问题；
它仍只允许已登记目标在 OpenMP 开启的配置中使用该选项。
基线与候选的原始诊断、分类及覆盖对照保存在本地 audit 目录，历史预算保持不变。

## 最终运行清单

| 配置/检查 | 结果 |
| --- | --- |
| GNU 15 Check，serial+MPI+JJ2LSJ OpenMP | 107/107 CTest |
| GNU 15 Release，同组程序；由 CMake 从显式 rCSFs 源目录构建 reader | 107/107 CTest |
| Release，MPI 关闭 | 69/69 CTest |
| Release，普通 serial apps 关闭、MPI 开启 | 85/85 CTest |
| rCSFs Rust 主库/原生 CLI/FFI及文档测试 | 217 passed，0 failed |
| 独立 native-reader | 4 passed，0 failed |
| 原 Python feature 编译 | 通过；未把旧安装扩展的测试冒充新源码的 Python wheel 验证 |
| 新增 Python 回归脚本 Ruff、Rust fmt、Git diff 检查 | 通过 |

RCI、GJ/RHFS/RIS4、biotransform、transition/phase 的串行小型产物与基线一致；
两个 rank 的对应 MPI 程序也完成基线产物对照。序列化记录和符号表检查确认原生 reader
仅导出所需 C ABI，运行库没有 Python/Arrow C++ 依赖。独立构建目录只随完整 Git 源码提供，
主 crate 归档包含共享实现和 C 头文件；CI 已补充原生 reader 测试步骤，未声称执行了远端 CI。

最终构建矩阵和测试数量以本地 `audit/native-validation.json` 及对应 CTest 日志为准。
该清单记录实际完成的配置，不把未运行的完整工具链组合或跨平台分发计入验收。

## V3 专项验收（2026-11-28）

Suite `5e5cfe6d9841666c7095f40995a439ecb40722f7` 在 macOS arm64 上重新构建，
**40/40 项专项 CTest 通过**，总用时 43.54 秒。构建采用 GNU Fortran 16.2.0、
`f2018-check` 和 Open MPI 5.0.11，启用串行应用、MPI 和 JJ2LSJ OpenMP。
依赖 rCSFs 的测试夹具由 4.1.0 原生 CLI 生成。[接入计划归档](../archive/plan-0006-csf-descriptor-v3-integration.md)

- 原生 Rust reader、C ABI 和 Fortran 适配覆盖块定位、损坏输入、失败后重开、描述符数组对照，以及计划指定的三条解析 CSF 和 bra/ket 交换。[原生 reader 集成测试](../../test/test_native_descriptors.py)
- 回归覆盖串行及 MPI 1、2、4 ranks，核对 MCP 结构和角系数、zero-first 块前缀、四种 RMCDHF 的 MCP/占据输入与标准权重，以及 RCI 的 V3 高 J 边界和续算。[MPI 描述符测试](../../test/rangular/test_mpi_descriptors.py)
- 双基组合并及端到端变换、跃迁、Landé 因子、JJ2LSJ 和 header-only 径向初猜均纳入本次选测。[CTest 登记](../../test/CMakeLists.txt)

从 Suite 根目录运行。`build_dir` 指向已完成构建的 Check 目录，配置时启用上述目标并通过
`RCSFS_TEST_CLI` 指定支持带索引 V3 的 CLI。

```fish
set -l build_dir ../temp/tests/suite-native-check
set -l test_root (realpath ../temp/tests)
set -lx TMPDIR (mktemp -d "$test_root/v3.XXXXXX")
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ctest --test-dir "$build_dir" --parallel 6 --output-on-failure \
    -R 'descriptor|state_workflow|single_subshell_tensor|rmcdhf_mcp_paths|rwfnestimate_header_input|^rci_restart|^lande_factors|^radiative_transitions|^jj2lsj_parallel_regression'
```

本机日志位于工作区的 `temp/tests/p6/build.log` 和 `temp/tests/p6/ctest-final.log`。
本次范围为本机专项回归。生产规模 MCP I/O 和资源预算继续按 [plan-0005](../plans/plan-0005-rmcdhf-mcp-production-io.md) 验收。
