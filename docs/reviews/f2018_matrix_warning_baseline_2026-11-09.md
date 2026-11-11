# F2018 构建配置矩阵与逐文件警告基线

后续状态：[优先警告修复第一批](f2018_warning_repairs_2026-11-09.md)已完成并更新基线，
当前为 32 配置、1,960 次 CTest 全通过及 475 条正确性复核诊断。
本文保留首次设施验收的 1,560 次测试和 515 条诊断统计；链接的基线文件现为更新版本。

日期：2026-11-09。源码基线：`f395146` 加此前
[strict 归一化时序修复（已移除）](../implemented/rejected_rmcdhf_variants.md)等未提交工作树。
本轮新增构建/诊断设施及回归；没有修改 Fortran 数值实现、输入、物理常数或数值门槛。

**状态：32/32 个配置通过，累计 1,560 次 CTest 执行全部通过，零失败/跳过；32 个配置的警告基线检查也通过。**
已登记的正确性警告不是已解决问题，不能据此宣布 P3/P9 或物理验收完成。

## 矩阵范围与执行契约

[matrix.py](../../test/f2018/matrix.py)枚举 Check/O0、Release/O2，分别取以下四个
Boolean 开关的全部 16 种组合，共 **32 个独立构建目录**：

| 标记 | CMake 选项 | 含义 |
| --- | --- | --- |
| s | `GRASP_BUILD_SERIAL_APPS` | 普通串行应用 |
| m | `GRASP_ENABLE_MPI` | MPI 库与应用 |
| c | `GRASP_BUILD_CSFG` | CSFG 专用库与程序 |
| o | `GRASP_JJ2LSJ_OPENMP` | jj2lsj OpenMP 后端 |

配置 ID 中 1/0 对应 ON/OFF。s=m=0 仍构建共享库、jj2lsj 和相应测试；
它不是“没有编译任何计算代码”的空配置。关闭 OpenMP 指 jj2lsj 后端，
不代表链接的数学库没有内部 OpenMP 支持。非线程专项固定 OMP/BLAS 为单线程。
MPI 专项按原测试覆盖 1/2/4 ranks；jj2lsj 线程专项保留原线程/逐位一致契约。

本次平台为 Linux x86_64、GNU Fortran 15.3.0。使用 Tools 共享 Python 环境；
`FC/OMPI_FC` 显式指向版本化 GNU 15，保持 Check/Release 的原严格选项。
两个配置并行、每个配置构建 `--jobs 3`；不是性能基准。
未复用旧 GCC 或其他配置的 `.mod/.o`，未安装或提交作业。

每个配置均核对实际 feature JSON、required 目标/测试、全部编译命令及 OpenMP 选项。
禁用功能时使用明确的相应子集，检查不应出现的目标/测试；不是将必需测试记为 skip。
矩阵使用默认容量、同一数学库及构建内准备工具；不涵盖任意用户 CMake 参数组合或外部准备工具。
`--case` 可复查子集，但不会得到 `full_matrix_acceptance=true`。

## 警告采集、基线与防回归

[warning_inventory.py](../../test/f2018/warning_inventory.py)作为 GNU Fortran 编译器 launcher，
为每次实际编译保存独立 JSON 记录，不从交错的并行构建文本中猜测警告归属。
只追加 JSON/无颜色诊断选项，保留原编译命令及已有 launcher。
每条记录包含命令、源文件/目标文件 hash 和原始结构化诊断；无警告编译也必须有记录。

严格 profiles 默认开启采集；[check.py](../../test/f2018/check.py)在运行 CTest 前强制核对：

- 实际编译源文件、目标集合/数量、编译命令、源码与 object hash 均有完整匹配记录。
- 基线编译器版本和配置匹配；新源文件/目标覆盖变化须显式审查。
- 警告身份包含目标、文件、行、列、选项和消息；各配置分别限制发生次数。
- 新诊断或次数增加失败，不能由另一个诊断/目标的减少抵消；`-w`、`-Wno-*` 不能隐藏警告。
- 缺失、陈旧、无法解析的记录失败；CTest 缺测试、失败、skip 等原 fail-closed 规则保持。

生成的[JSON 基线](../reference/f2018/gfortran15-warning-baseline.json)和
[逐文件索引](../reference/f2018/gfortran15-warning-baseline.md)只保存仓库相对路径、
源码 hash 与诊断，不包含构建目录或机器路径。完整原始记录保留在忽略的构建目录。
JSON 基线覆盖 **3,040 个源文件、4,841 条唯一诊断**；按目标区分，同一源位置在不同
目标中的相同消息算不同诊断。923 个文件有警告，其余文件也登记在覆盖清单中。

| 分类 | 唯一诊断数 | 含义 |
| --- | ---: | --- |
| correctness-review | 515 | 未初始化、越界、截断、类型转换、接口漂移、隐式静态数组等候选风险 |
| obsolescent-standard | 403 | 编译器标记为过时但仍属标准的语法，不等同于非法扩展 |
| general | 3,923 | 其余诊断；并不意味着已经逐条判断为无害 |

**195 个文件包含正确性复核项。** 例如普通 RCI 串行/MPI 的 `NCTEI` 未初始化诊断、
Davidson 手写接口 INTENT 不一致、CSFG DACON 的固定数组下标、字符截断及大局部数组
自动转为静态存储。它们需后续逐项判定并建立触发/oracle；本轮没有把警告逐条裁定为真缺陷。
全功能 Check/Release 分别记录 4,621/4,873 次警告，正确性复核项为 340/504 次。
Release 的优化分析暴露更多可能未初始化诊断，不能以 Check 清单替代 Release。

`warning_gate.status=passed` 只表示没有超过固定基线；
`correctness_warnings_clear=false` 仍明确保留。计划中相关诊断必须清零的 P3/P9 条件没有放宽。
基线不会自动增长；显式生成另一份候选后审查差异，包括修复后的减少，再有意更新基线。
链接器警告与测试脚本内部临时编译探针不在此 CMake Fortran 编译基线内，仍保留于相应日志。

## 验证设施本身的检查

新增 `f2018_matrix_warnings_driver`（7 个 Python 检查），覆盖 32 个组合、解析失败、
编译器版本变化、零警告源文件遗漏、目标覆盖变化、陈旧源码/object/命令、schema 变化、
新警告不能抵消以及单条次数增加。原 acceptance 驱动扩展了 OpenMP 实际命令检查和 `-w` 反例。
真实 GNU 15 launcher 探针验证带空格目录、嵌套构建路径和已有 launcher 的保留；
对真实采集结果将一个诊断预算减一，确认命令返回失败并报告新增/超额诊断。

首次采集中编写 launcher 时有一次并发读取尚未写完文件的失败；代码固定后保留失败日志，
重跑受影响配置并重新核对全部 32 个构建。正式基线来自完整记录，未把这次失败当作 Fortran 缺陷或通过。

## 矩阵验收结果

两种 profile 各 16 个配置、780 次 CTest 执行，全部通过；累计核对 56,624 条 Fortran 编译记录。
每个配置的 `strict_acceptance=true`，总表 `full_matrix_acceptance=true`；同时每个配置仍有
`correctness_warnings_clear=false`。这里的通过范围是当前目标/测试和固定警告预算，
不等于全仓库已无未定义行为或完整迁移验收。

| 配置（s/m/c/o） | 编译单元 | CTest 通过 | 警告次数 | 正确性复核次数 |
| --- | ---: | ---: | ---: | ---: |
| `check-s0-m0-c0-o0` | 528 | 29/29 | 345 | 93 |
| `check-s0-m0-c0-o1` | 528 | 29/29 | 344 | 93 |
| `check-s0-m0-c1-o0` | 962 | 34/34 | 975 | 185 |
| `check-s0-m0-c1-o1` | 962 | 34/34 | 974 | 185 |
| `check-s0-m1-c0-o0` | 1564 | 56/56 | 1373 | 154 |
| `check-s0-m1-c0-o1` | 1564 | 56/56 | 1372 | 154 |
| `check-s0-m1-c1-o0` | 2484 | 67/67 | 3818 | 292 |
| `check-s0-m1-c1-o1` | 2484 | 67/67 | 3817 | 292 |
| `check-s1-m0-c0-o0` | 1298 | 30/30 | 1149 | 141 |
| `check-s1-m0-c0-o1` | 1298 | 30/30 | 1148 | 141 |
| `check-s1-m0-c1-o0` | 1732 | 35/35 | 1779 | 233 |
| `check-s1-m0-c1-o1` | 1732 | 35/35 | 1778 | 233 |
| `check-s1-m1-c0-o0` | 2334 | 64/64 | 2177 | 202 |
| `check-s1-m1-c0-o1` | 2334 | 64/64 | 2176 | 202 |
| `check-s1-m1-c1-o0` | 3254 | 75/75 | 4622 | 340 |
| `check-s1-m1-c1-o1` | 3254 | 75/75 | 4621 | 340 |
| `release-s0-m0-c0-o0` | 528 | 29/29 | 399 | 107 |
| `release-s0-m0-c0-o1` | 528 | 29/29 | 398 | 107 |
| `release-s0-m0-c1-o0` | 962 | 34/34 | 1090 | 216 |
| `release-s0-m0-c1-o1` | 962 | 34/34 | 1089 | 216 |
| `release-s0-m1-c0-o0` | 1564 | 56/56 | 1485 | 224 |
| `release-s0-m1-c0-o1` | 1564 | 56/56 | 1484 | 224 |
| `release-s0-m1-c1-o0` | 2484 | 67/67 | 4027 | 413 |
| `release-s0-m1-c1-o1` | 2484 | 67/67 | 4026 | 413 |
| `release-s1-m0-c0-o0` | 1298 | 30/30 | 1246 | 198 |
| `release-s1-m0-c0-o1` | 1298 | 30/30 | 1245 | 198 |
| `release-s1-m0-c1-o0` | 1732 | 35/35 | 1937 | 307 |
| `release-s1-m0-c1-o1` | 1732 | 35/35 | 1936 | 307 |
| `release-s1-m1-c0-o0` | 2334 | 64/64 | 2332 | 315 |
| `release-s1-m1-c0-o1` | 2334 | 64/64 | 2331 | 315 |
| `release-s1-m1-c1-o0` | 3254 | 75/75 | 4874 | 504 |
| `release-s1-m1-c1-o1` | 3254 | 75/75 | 4873 | 504 |

## 复现与证据

常规验证使用已登记基线：

```sh
export FC=/usr/bin/gfortran-15 OMPI_FC=/usr/bin/gfortran-15
../graspkit-tools/.venv/bin/python test/f2018/matrix.py build-gf15-matrix-20261109 \
  --parallel-configs 2 --jobs 3
```

初次采集/显式候选流程见[测试说明](../../test/f2018/README.md#complete-configuration-matrix-and-warning-baseline)。
旧的无记录构建须使用新目录或 `cmake --build <build> --clean-first`；
仅有增量构建文本日志不能证明无警告文件或旧 object 已覆盖。

本地证据位于忽略的 `build-gf15-matrix-20261109/`：

- `matrix-build.json`、`matrix-validate.json`：完整矩阵及逐配置结果。
- `<case>/audit/acceptance/`：配置、实际编译命令核对、构建、JUnit、当前源码 hash。
- `<case>/audit/warning-records/`、`warnings.json`：逐 object 原始记录、归并诊断。
- `<case>/audit/matrix-case.json`：测试/告警结果及验收摘要 hash。
- `launcher-probe.json`、`warning-budget-negative.log`：真实 launcher 与拒绝超额诊断的反例。
- `matrix-build-first-attempt.json` 及失败配置的 `audit/launcher-edit-failure.log`：首次采集失败证据。

没有修改科学计算参数、放宽数值容差、install、提交、push 或 SLURM 作业。
完整 F00–F10、真实多根/属性/跃迁、接口/状态清理、mpi_f08/大计数及物理验证仍按计划推进。
