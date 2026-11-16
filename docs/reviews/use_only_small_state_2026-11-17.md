# USE ONLY 第四批：小型共享状态模块

日期：2026-11-17。基准为 HEAD `3e02e0f` 加此前未提交修改，接续
[共享库过程导入整理](use_only_library_interfaces_2026-11-17.md)。本批保持接口、精度和计算行为。

## 修改范围

220 份源码补齐 284 条 ONLY。本批选择不依赖其他模块、公开变量数量不超过三个的
13 个状态模块，并用修改前 GNU 编译模块的导出符号表确认完整列表。
ONLY 列出原有全部公开变量，不依据局部使用次数删减，不修改变量声明或初始化值。

| 模块 | ONLY 列表 | 修改数量 |
| --- | --- | ---: |
| iounit_c | stderr_unit, stdin_unit, stdout_unit | 153 |
| foparm_c | iccut | 37 |
| orthct_c | orthst | 20 |
| iccu_c | iccut, iccutblk | 20 |
| offd_c | noffd1, noffd2 | 9 |
| orbord_c | nordff, nordii | 8 |
| core_c | closed_core_subshell_count | 8 |
| kkstartbreit_c | kmax, kstartbreit1, kstartbreit2 | 8 |
| jlabl_c | jlbl, jlbp, jlbr | 7 |
| fposition_c | nposition | 6 |
| qedcut_c | nqedcut, nqedmax | 4 |
| mtjj_c | mt | 2 |
| mtjj2_c | mt11, mt9 | 2 |

范围为 `src/` 内外部程序单元或单模块文件 CONTAINS 后的过程内 USE，
未改变模块头部重导出、带重命名的导入、测试预处理文件或复杂状态模块。
逐行核对只在目标 USE 后追加 ONLY，行数不变，没有新增超长行。

94 份导入 RMCDHF 文件追加
[状态模块 ONLY 来源快照](../reference/rmcdhf_use_only_state_snapshot.json)，绑定之前的
RMCDHF USE ONLY manifest；历史快照不覆盖，来源检查继续验证完整前驱链和精确哈希。

## 验证

使用本机 Homebrew GNU Fortran 15.3.0，显式设置 FC/OMPI_FC；服务器工具链约定不变。
Python 使用工作区共享环境，采用独立 Release/Check 构建目录。

- 13 个模块的编译器导出表均与完整 ONLY 列表一致，未减少任何公开变量。
- Fortitude C121 从 3,063 降至 2,779，默认诊断从 7,119 降至 6,835。
- 扩展正确性检查从 14,551 降至 14,267，差异恰为 284 条 C121；其余诊断按文件、规则、
  行号和消息比较完全一致。
- 220 份文件的编辑范围和修改后哈希核对通过；源码与历史 manifest 的来源链篡改负例通过。
- Release/Check 全量构建通过。256 个受影响 Release 编译对象逐字节一致，全部 1,204 份
  编译模块文件解压内容一致，接口未变。所有修改源码均列入当前构建目标，重复编译分别计数。
- 492 份导入 RMCDHF 文件的完整来源链及精确哈希检查通过。
- Release/Check 两套完整 CTest 均一次性 89/89 通过；`git diff --check` 通过。

本地忽略目录 `build-use-only-state-20261117/` 保存修改前源码、每条导入列表、预期与编译器
实际导出表、前后哈希、静态诊断以及编译对象/模块接口对照和回归日志。
没有新增物理计算；既有科学验证范围和 macOS memory 完整 SCF 限制仍见
[S3/S4 报告](static_numeric_interfaces_2026-11-17.md)。未提交、推送或安装覆盖生产二进制。

剩余 2,779 条 C121 继续按符号来源、别名、泛型和重导出契约分批核定。
