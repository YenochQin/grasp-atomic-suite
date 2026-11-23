# S1–S4 后 fortls / Fortitude 全量复查

日期：2026-11-17。基准 HEAD：`3e02e0f407d085504d9b2fb8445d50e8d013b37c`；
检测对象包含当时尚未提交的 S1–S4 源码与测试，后续纳入 `bea6655`。

> 本文保留工具解析/编译器对照证据。后续修复及 CSF 输入迁移已实施，
> 当前待办见[剩余事项](../plans/plan-0004-f2018-follow-up-validation.md)，下文结论仅对应扫描基准。

结论：本次未发现这些修改新增的已确认功能缺陷，但仍有属于描述符迁移任务的 CSF 输入风险，
以及大量风格、接口规范和工具解析覆盖问题。本次只检测与记录，没有修改生产源码、自动修复、提交或推送。

## 工具与覆盖

- 直接使用用户通过 `uv tool install` 安装的 **fortls 3.2.2**、**Fortitude 0.10.1**。
  没有安装或升级工具，没有把工具依赖装入共享运行环境；fortls 自动更新关闭。
- 当前输入为 **1,731** 份文本：`src/` 1,697 份（含 `memory_man.f90_PROFILING`）、
  `test/` 34 份。按显式文件清单扫描，不包括 build 副本和 vendor 的 lint 诊断。
  生成在 Python 字符串中的 Fortran 驱动不算作独立源文件。
- Fortitude 分别运行默认规则和 `C,PORT,OB,E` 加 preview，使用 Fortran 2018、132 列，
  `--isolated --no-fix --no-unsafe-fixes`。未运行全规则 ALL 扫描，不宣称 S252 等旧工具限制已解决。
- fortls 使用实际 CMake 编译命令、宏和目标依赖：**78 个编译目标上下文 + 22 个未构建/动态文件补充上下文**。
  通过已安装工具的 `LangServer.file_init`、链接解析和 `get_diagnostics` 引擎执行。
  不把不同程序的同名模块放入一个全仓库符号表；100 个上下文均无加载/诊断引擎异常。
- 当前 HEAD 另导出到忽略目录，用同工具、同配置重新扫描，作为本次差异对照。
  HEAD 有 1,727 份对应文本，当前多四份 Fortran 测试文件。不是拿旧日期的诊断总数直接比较。
- 扫描前后 1,731 份源码 SHA-256 全部一致。

## Fortitude 结果

| 口径 | 同工具扫描 HEAD | 当前工作树 |
| --- | ---: | ---: |
| 默认规则，132 列 | 30,993 | 31,061 |
| Correctness / Portability / Obsolescent / Error，加 preview | 16,351 | 16,377 |
| C181：IOSTAT 检查提示 | 13 | 10 |
| E001：解析语法提示 | 19 | 20 |

默认规则当前分类：Style 20,804、Correctness 8,024、Modernisation 1,316、
Obsolescent 897、解析 Error 20。规则类别或工具标注的 error 不等于编译错误，也不等于独立功能缺陷。

IOSTAT 的三个消失位置分别是 `ris4/angdata`、`jj2lsj::uniquelsj` 和 MPI worker 输出重定向。
当前 C001/C002/C021/C022/C101/C161 均为零；同一 HEAD 基线也已为零，
不能将其全部归功于本批修复。尤其显式 `_real32` 不再触发缺少 kind 的提示，
因此 C021 为零不能替代此前的双精度边界验证。

### 剩余 10 条 C181

| 位置 | 数量 | 当前核定 |
| --- | ---: | --- |
| `lib9290/lodcsh.f90:97`、`lodcsl.f90:177` | 2 | 已知 CSF 输入风险，属于描述符迁移任务，未因本批工作关闭 |
| 三份 `rtransition*/lodcslm.f90:163` | 3 | 已知 CSF 记录状态被覆盖风险，同上 |
| `rangular/rangular_input.f90:76` | 1 | BLOCK 外紧接着检查状态，工具的同作用域限制 |
| `npol_init.f90:99,112,119` 与 `npol_input.f90:49` | 4 | 清理错误状态不覆盖原失败，或后续 READ 已受状态保护；保留既有复核分类，不称为错误输入被接受 |

NP 的 INQUIRE 输出引用等仍可作为独立防御性改善；本次没有把静态提示数量当作实际缺陷数。

### 新增提示与解析覆盖

与 HEAD 比较，新增提示主要是：

- 检查 IOSTAT 后尚未提供 IOMSG 的 C183；生产代码仍通过既有处理器报告阶段、unit 和状态码。
  C183 总数净增 14，不等于新增 14 个未检查错误。
- JJ2LSJ 两条跨行单行 IF 的 C152 可读性建议。
- 新测试的 USE/ONLY、固定 unit、数组长度常量和 OPEN action 建议。
- `test/ris4/cache_io.F90` 的条件预处理 ELSE 被误解析，新增一个 E001。
  GNU 15 按实际串行/MPI 宏预处理后，两份输出分别重新运行 Fortitude E001，均为零。

其余 E001 位于既有续行 FORMAT、算术表达式和条件预处理源码。
对 E001 涉及文件及 fortls Error 所在文件，按当前实际编译条件重新运行 GNU 15：
**21 份源码、28 个编译上下文全部通过**。
这排除了这些报告作为当前编译语法错误的证据，但 Fortitude 对原文本的 AST 覆盖仍可能不完整。

## fortls 结果与对照

原始扫描得到 53 条独立诊断。为动态编译的 `interface_capacity.f90` 补齐
其测试脚本明确列出的 `tail_I/frmtfp_I/frmrwf_I` 三份依赖后，最终为：

| 类别 | HEAD | 当前 |
| --- | ---: | ---: |
| Error | 3 | 3 |
| Warning | 28 | 28 |
| Information | 18 | 19 |
| 合计 | 49 | 50 |

按文件、严重性和消息比较（忽略改动造成的行号移动），唯一新增项是新缓存 MPI 测试中的
`Module "mpi" not found in project` Information。

- **3 条 Error**：均在 `csf_descriptor_input.f90`，仍是 `block_lengths(...)` /
  `block_labels(...)` 数组赋值被误认为 BLOCK，导致后续过程的 Invalid parent。
  本次用实际 fortls CLI 重新运行最小对照：`block_labels(1)=1` 产生 Invalid parent；
  仅把临时样本变量改成 `labels` 后消失；两份样本均通过 GNU 15。
  没有为规避解析器而重命名生产变量。
- **28 条 Warning**：27 条宿主名称遮蔽（其中 7 条由上述伪嵌套带出），
  另有 `ris_cal.f90` 一条超过 132 列的提示；均未新增。
- **19 条 Information**：全部是源码项目中找不到外部 MPI 模块。
  MPI 模块由已验证的编译工具链提供；不能据此判定 MPI 安装/链接失败。

旧报告中两份 `vpint` 的错误未在本次出现；这里只描述当前扫描结果，
不把它解释为整个 fortls 解析器的 fixed/free-form 问题已经解决。

## 证据与后续工作

本地忽略目录 `build-static-rescan-20261117/` 保存：

- `scan-files.json`、`source-hashes.json`：当前完整清单及源码身份；
- `current-{default,correctness}.json`、`head-{default,correctness}.json`、对应 stderr 和摘要；
- `fortls-contexts.json`、`fortls-results.json`、`fortls-final-unique.json` 及 HEAD 对照；
- `fortls-dynamic-control.json`、`parser-probes/`：依赖补齐和 CLI 误报对照；
- `compiler-controls.json`、`cache-*-preprocessed.f90`、`cache-*-lint.json`：编译器和预处理对照；
- `scan_fortitude.py`、`scan_fortls.py`、`fortls_controls.py`、`compiler_controls.py`：本次执行脚本。

当时建议通过描述符迁移处理五处 CSF 输入风险；该迁移现已完成，见[原生输入说明](../implemented/design-0001-native-csf-descriptors.md)。
风格、显式接口规范、无可达缺陷证据的可维护性建议另行分批；不以全库自动修复换取零诊断。
本轮未重新运行 SCF 或全量 CTest；已有软件/数值验证边界见[兼容性汇总](../implemented/design-0004-f2018-compatibility.md)，工具扫描不能替代这些验证。
