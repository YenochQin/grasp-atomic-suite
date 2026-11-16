# 非 CSFG Fortran 源码小写化与验证记录

日期：2026-11-16。基准：`caaf3a51ce12ec8e85e30622fc854c9e2a070683`，叠加本次工作树改动。

状态：普通 Fortran 源码已完成大小写统一；Release 和 Check 均构建成功，软件回归通过。
没有实施语义重命名，没有修改数值表达式、接口含义、数组布局或文件格式。
这次验证不构成完整原子计算的物理精度或收敛性验收。

## 范围与转换规则

采用[变量清单](../plans/non_csfg_variable_inventory_2026-11-16.md)的范围：

- `src/` 下 1,879 个 `.f90` 文件及 `src/lib/libmod/memory_man.f90_PROFILING`。
- 排除路径中含 `csfg` 的文件，以及仅由 [CSFG 构建入口](../../cmake/CSFGExtension.cmake)
  启用的 `src/tool/graspg/`、`src/lib/libdvd90mpi/`。
- 不改文件名；不转换 CMake、Makefile、Shell 文件；`test/` 下的 Fortran 测试源码也未批量转换。
- Fortran 关键字、标识符、数值指数标记及未加引号的 FORMAT 描述符转换为 ASCII 小写。
  字符串（含续行及成对引号）、注释、OpenMP 指令行、预处理行和宏引用保留原字节。
  空格、缩进、换行、标点和非 ASCII 字节不变。

使用临时词法转换器先完整预览、保存原文，再写入。宏保护名单来自 GNU Fortran 15 的
预定义宏、编译命令的 `-D` 参数和源内预处理声明，共 49 个名称。
当前纳入的源码没有预处理行或这些宏的实际引用；宏、预处理续行、引号、字符串续行、
插入注释行、CRLF 和 BOZ 字面量的保护行为另用一次性小样本验证。
遇到未处理的 Hollerith 或未闭合字符串会拒绝转换；本次没有此类阻断项。

| 项目 | 结果 |
| --- | ---: |
| 纳入 Fortran 文本 | 1,880 |
| 实际修改文件 | 1,872 |
| 原已符合要求的文件 | 8 |
| 改变大小写的词法片段 | 269,251 |
| 大写变小写的字节数 | 963,036 |
| 原样保留的字符串片段 | 13,901 |
| 原样保留的整行注释 | 50,924 |
| 原样保留的行内注释 | 438 |
| `src/` 中范围外、路径和哈希均不变的文件 | 1,006 |

逐字节验证只允许 `A`–`Z` 变成对应小写字母，并核对全部受保护片段；文件长度不变，
转换幂等。fparser 0.2.5 对 1,880 对前后文本解析成功，忽略代码大小写后的语法树表示一致；
字符串和注释的精确一致性由独立字节检查保证。`storage_bytes.f90` 的 Fortran 2018 assumed-rank
仅在解析器内存副本中作适配，源码仍保留原来的 `(..)`。

## 编译与回归

工具链为 GNU Fortran 15.3.0、CMake 4.5.0、系统 Open MPI/OpenBLAS；显式设置
`FC=/usr/bin/gfortran-15`、`OMPI_FC=/usr/bin/gfortran-15`，MPI wrapper 返回同一编译器。
Python 使用工作区共享环境。构建使用 Unix Makefiles；没有安装或覆盖现有 `bin/`。

全部配置均启用普通串行、MPI、jj2lsj OpenMP 和 rangular OpenMP，CSFG 关闭。
基线先从未转换的工作树完成编译和测试，并保存二进制；转换后使用另外两个全新构建目录。

| 配置 | 编译 | CTest |
| --- | --- | --- |
| 转换前 `f2018-release` | 成功 | 67/67，通过，无跳过 |
| 转换后 `f2018-release` | 成功 | 首轮 64/67；三个维护检查适配后定向重跑 3/3，最终 67 项均通过 |
| 转换后 `f2018-check` | 成功 | 67/67，通过，无跳过 |

Check 启用既有的 Fortran 2018、显式接口、运行时边界、回溯与浮点异常检查选项。
回归包括共享网格、Davidson、CI 阻尼、径向积分、核极化组件、描述符输入、MPI 通信、
OpenMP 角系数、密度/自然轨道、jj–LS 转换和高角动量轨道初值等已有测试。
未注册到普通 CMake 的旁存源码由词法/语法树检查覆盖，不能称为逐文件独立编译验收。

首轮失败项均与源码文本形式有关，数值回归没有执行失败：

1. [网格 patcher 测试](../../test/test_patch_grasp_grid.py)的预期替换及负例构造原来假定大写，
   现按 Fortran 不区分大小写匹配；仍精确比较其余字节及备份恢复结果，并检查负例确实删去一个重导出名。
2. [容量接口回归](../../test/standard_grasp/alternate_capacity.py)改为匹配完整的参数声明，
   保留原拼写，要求两个默认容量各恰好替换一次；没有放宽容量值断言。
3. [上游来源检查](../../test/npol/check_upstream.py)加入独立的
   [小写化快照](../reference/rmcdhf_lowercase_snapshot.json)。480 份 Fortran 导入文件的转换前原文，
   已用原检查逻辑逆转已审查的 F2018/NP 改动并核对上游 SHA-256；随后记录转换前、转换后的精确 SHA-256。
   新快照同时固定原上游 manifest 和 F2018 compatibility manifest 的哈希。
   **原始两个来源 manifest 没有修改**，其他 12 份导入文件沿用原检查。
   这里不是把全文小写后再比较哈希：当前 Fortran 文件必须逐字节等于已审查的小写化版本。
   在一次性副本中篡改注释、字符串或代码的三个负例均被拒绝。

## 目标文件与输出对照

转换前后 Release 均产生 2,227 个目标文件；其中 2,013 个完全相同。
其余 214 个的差异只位于 `.rodata.str1.1` / `.rodata.str1.8`，合计 2,793 个字节，
全部是 ASCII 大写转小写；所有其他 ELF section 均一致，包括指令、重定位及数值数据。
原因是未加引号的 FORMAT 语法被编译器存入只读字符串，例如
`(1P,6(1X,1D19.12))` 变为 `(1p,6(1x,1d19.12))`；源码中原本加引号的文本未变。
25 个应用二进制中只有 1 个整体哈希相同，因此**不宣称全部可执行文件字节一致**。

又使用现有测试的合成输入直接运行转换前/后的二进制，不采用新的物理模型或精度阈值：

| 用例 | 条件 | 精确一致的输出文件对 |
| --- | --- | ---: |
| rangular | 串行、同一 CSF 输入 | 6 |
| rangular_omp | 4 线程、小分块缓存 | 6 |
| rangular_mpi | 2 ranks、同一 CSF 输入 | 12 |
| rdensity | 加权多块及有主量子数间隔两组输入 | 4 |
| rdensity_mpi | 同上，2 ranks | 4 |
| jj2lsj | even、odd、f3 三组，4 线程，完整转换 | 12 |
| 合计 | 10 组用例 | 44 |

MCP、密度/自然轨道文件、带符号转换系数和标签等输出全部逐字节一致。
最初一次附加 MPI 对照使用过长工作路径，未转换的基线程序在 `sys_getwd` 即停止；
改用短临时路径后两版均完成。该问题属于既有路径长度限制，没有因此修改生产代码。

## 保留的限制

- 历史 GNU 15 告警预算未通过源文件覆盖检查：本次 Check 编译 1,911 个不同源码文件，
  预算覆盖 1,874 个，缺少 37 个。转换前 Release 的实际编译源码清单与转换后完全相同，
  证实此覆盖差异先于本次小写化。没有更新或放宽历史预算。
  Check 采集到 2,024 条告警，其中 243 条被现有分类器归为 correctness-review；
  不能据此宣称告警清零、完整告警门禁通过或这些风险已解决。
- `git diff --check` 指出 10 处既有行尾空白。与转换前对应行核对，空白完全一致；
  为保持纯大小写改动，没有顺便清理这些行。
- 未运行 CSFG、完整 36 配置矩阵、生产规模 SCF 或物理收敛实验。
  上述结果支持这次大小写变换的软件行为保持，不扩大为所有算法或物理结果的验证。

## 本地追溯文件

所有构建、原文备份、一次性转换/验证脚本及数值输出均保存在被 Git 忽略的
`build-lowercase-20261116/`，没有提交计算输出或本机环境：

- `source-before/`、`source-before.json`：转换前原文及全部 src 哈希。
- `lowercase-{preview,apply,verify}.json`、`lowercase_fortran.py`：保护规则及逐文件变化。
- `ast-equivalence.json`、`validate_case_ast.py`：1,880 对语法树比较。
- `release-baseline/`、`release-after/`、`check-after/`：三套独立构建。
- `logs/ctest-baseline.log`、`logs/ctest-release-initial.log`、`logs/ctest-release-rerun.log`、
  `logs/ctest-check.log`：完整测试日志；各构建目录保留对应 JUnit XML。
- `artifact-equivalence.json`、`object-section-comparison.json`：目标文件及二进制比较。
- `differential-summary.json`、`differential-verified/`：44 对精确输出及运行日志。
- `warning-baseline-gap.json`、`check-after/audit/warnings.json`：历史预算覆盖差异与实际告警。
- `whitespace-audit.json`、`final-audit.json`：原有空白核对及最终检查汇总。

后续语义重命名继续依据变量清单逐作用域核实，不能把这次 `NCFTOT → ncftot` 当作
已经完成 `ncftot → total_csf_count`。
