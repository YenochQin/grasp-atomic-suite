# 风格第三批：小数写法

日期：2026-11-17。基准为 HEAD `3e02e0f` 加未提交的 S1–S4、空格及
[引号/END 整理](style_quotes_end_2026-11-17.md)。本批仅处理 Fortitude S291。

## 修改范围与精度约束

153 份 Fortran 文件中的 656 条提示已处理，S291 从 672 条降至 16 条。
例如 `1.e9_real64` → `1.0e9_real64`、`.998_real64` → `0.998_real64`。
只给小数点缺失的整数部分或小数部分补零，不改指数、符号、kind、变量类型或表达式。
特别是原有 `_real32` 不改为 `_real64`，避免把风格整理变成精度变更。

只应用 S291 标记 safe 的插入零操作；逐处验证十进制值、指数原文、kind 原文一致，
除此之外源码字节保持不变。物理行数不变，没有新增超过 132 列的行。
剩余 16 条均在四份报解析错误的 `strsum.f90` 中，整份文件保留。
未运行全仓库自动修复或使用 unsafe fixes。

76 份导入 RMCDHF 文件追加
[decimal 来源快照](../reference/rmcdhf_decimal_snapshot.json)，绑定上一份 quote/end manifest；
历史 manifest 不变，来源检查继续验证前驱及每份文件的精确哈希。

## 验证与证据

使用本机 Homebrew GNU Fortran 15.3.0，显式指定 FC/OMPI_FC；未修改仓库的服务器
工具链路径约定。Python 使用工作区共享环境，采用独立 Release/Check 构建目录。

- 所有修改覆盖 49 对不同常量写法。GNU Fortran 编译运行探针，逐对核对 KIND 与
  TRANSFER 到字节数组后的位模式，全部一致；不使用浮点误差容差。
- Release/Check 全量构建通过；对照上一批 Release 构建，受影响源码对应的
  178 个编译对象全部逐字节一致。多目标重复编译分别计数，未编译文本不计入对象对照。
- 两套 CTest 首轮各 88/89 通过；唯一失败是数值边界测试的源码提取器只匹配
  `.998_real64`。提取器现接受可选前导零，原有阈值、出现次数和邻接浮点数边界断言
  不变。两套定向复测各 1/1 通过，合计各 89 项通过。
- Fortitude 默认诊断从 11,567 降至 10,911，恰好减少 656 条。
- correctness/portability/obsolescent/error 加 preview 的 16,377 条诊断，按文件、规则、
  行号和消息比较完全一致。
- 492 份导入文件来源检查通过；来源链篡改负例通过，包含修改源码及 declaration、numeric-kind、spacing、quote/end
  历史 manifest 的情况。

本地忽略目录 `build-style-decimal-20261117/` 保留修改前源码、逐处编辑记录、前后哈希、
限定规则执行脚本、常量位模式探针、静态诊断和构建/测试日志。
本批不开展新的物理计算；已有科学验证范围与 macOS memory 完整 SCF 限制仍见
[S3/S4 报告](static_numeric_interfaces_2026-11-17.md)。未提交、推送或安装覆盖生产二进制。
