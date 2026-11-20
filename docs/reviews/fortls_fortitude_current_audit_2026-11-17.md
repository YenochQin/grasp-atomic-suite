# Fortitude / fortls 当前 HEAD 审查

日期：2026-11-17。基准：`5806884a6c448fab7c6bc8ca4ccfe84908bbd5b8`，开始时工作树干净。
本轮重新执行静态工具、复核重点源码并记录结果；未修改生产源码、自动修复、提交或运行原子计算。

> 2026-11-17 文档整理注：本文的“当前”指扫描基准 `5806884`。随后 `8cf6b1d`、
> `80f6e76` 已替换生产 CSF 文本输入；下述五处旧 reader 风险不再列为当前生产待办。
> 原始扫描结果与限制保留，迁移后的诊断需重新扫描，见[接续清单](../plans/static_diagnostics_remaining_2026-11-17.md)。

## 需要优先修复的输入风险

1. **P1：首条 CSF 读取失败后继续解析。** `src/lib/lib9290/lodcsh.f90:97` 的 READ 设置
   `ios` 后立即调用 `prsrcn`，没有检查读取状态；其返回的 `ierr` 也未检查就计算电子数。
   缺失或损坏首条 CSF 时，后续计数没有可靠输入保证。应先验证读取和解析状态，再累加占据数。
2. **P1：检查读取状态之前使用记录，可能覆盖原始错误。**
   `src/lib/lib9290/lodcsl.f90:178` 和
   `src/apps/rtransition/lodcslm.f90:164`、`rtransition_mpi/lodcslm.f90:164`、
   `rtransition_phase/lodcslm.f90:164` 都先检查 `record(1:2)`，遇到块分隔符时再 READ 并覆盖 `ios`，
   最后才检查状态。异常/EOF 路径不应使用失败读取后的记录，也不应把原始错误交由第二次读取覆盖。
   应分别处理正常记录、允许的 EOF 和读取错误，再识别分隔符。

这五处是现存风险，属于[CSF 描述符迁移](../implemented/native_csf_descriptors.md)。
本轮依据当前源码确认控制流问题，没有重新注入损坏文件，不宣称已复现错误物理结果。

## 工具结果

工具版本：Fortitude 0.10.1、fortls 3.2.2。覆盖 src/test 共 **1,731** 份 Fortran 文本，
含 `memory_man.f90_PROFILING`；vendor 只作为 fortls 依赖，不计入第一方诊断。
扫描前后全部输入 SHA-256 一致。

| 检查 | 结果 |
| --- | ---: |
| Fortitude 默认规则，F2018、132 列 | 4,057 |
| Fortitude C/PORT/OB/E，加 preview | 11,489 |
| fortls Error | 3 |
| fortls Warning | 28 |
| fortls Information | 18 |

两个 Fortitude 口径不能相加；均使用 isolated、no-fix、no-unsafe-fixes。
默认规则的主要项是 C061 缺少 INTENT 1,104、C131 可见性策略 803、OB031 783、
C092 677。C121 只剩 `src/lib/mpi90/mpi_C.f90:6` 一条，是 MPI 模块导入；
不能为清零提示随意收窄重导出接口。C001/C002/C021/C022/C101/C161 本轮均为零。
规则提示数不是独立缺陷数，未执行 ALL 全规则扫描。

C181 共 10 条：上述五处为输入风险；rangular_input 的状态在 BLOCK 外紧接着检查；
NP 四处涉及保留原始错误的清理、受保护的后续读取或 INQUIRE 后的防御性检查，
不把这五条剩余提示判为已证实的错误输入接受。

## fortls 上下文与解析限制

新建 GNU Fortran 15.3 / MPI 的 CMake 配置，启用串行应用并导出编译命令；没有执行完整构建。
按生成的目标依赖、include 和宏建立 **78 个目标上下文、22 个补充上下文**，
避免同名应用模块混入一个符号表。通过已安装 fortls 的 LangServer 解析/链接/诊断引擎执行，
禁用自动更新；100 个上下文无引擎异常，1,731 份输入全部检查。
原始 52 条诊断，按动态测试脚本补齐 interface_capacity 的依赖后为 49 条。

- 3 条 Error 均位于 csf_descriptor_input 的后续过程声明，符合既有 block_* 数组元素赋值被
  误认作 BLOCK 的解析问题。本轮 CLI 最小对照再次复现：`block_labels(1)` 赋值引出
  Invalid parent，改成 `labels(1)` 后消失，两份样本均通过 GNU 15 F2018 语法检查。
- 28 条 Warning 为 27 条宿主名称遮蔽及一条行宽提示；部分遮蔽随伪嵌套产生。
- 18 条 Information 全部为源码项目中找不到外部 mpi 模块。CMake 本轮成功发现 MPI；
  这些提示本身不是 MPI 链接失败的证据。

Fortitude 有 20 条 E001，位置仍是续行 FORMAT、表达式和带预处理条件的测试。
这些模式在[前次复查](fortls_fortitude_rescan_2026-11-17.md)有编译器核查记录；
**本轮未对当前全部 E001 文件重做编译器验证**，不将历史验证当作当前编译通过，
也不直接把 E001 当作源码语法缺陷。受影响文件的 AST 检查覆盖存在限制。

## 证据和边界

本地忽略目录 `build-static-audit-20261117/` 保存本轮脚本、输入清单与哈希、
两份 Fortitude JSON、fortls 上下文/结果/最终去重诊断、动态依赖补齐结果、CLI 解析最小对照和日志。
`cmake/` 保存本轮新生成的编译上下文。原始输出和本机路径不纳入提交。

优先完成五处 CSF reader 的错误处理或迁移，再分批核定 INTENT 和公共模块可见性。
本轮未运行全量 CTest、SCF、MPI 数值回归或物理收敛验收；静态工具不能替代这些验证。
