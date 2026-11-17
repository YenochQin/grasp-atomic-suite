# USE ONLY 第九批：数组状态与常量

日期：2026-11-17。基准为 `bea6655` 加未提交的第七、八批修改。本批未提交、未推送。

146 份文件补齐 193 条 ONLY，覆盖 prnt_c、cons_c、mcpdata_c、m_c、bess_c、syma_c、
trk_c、int_c 八个模块。完整保留原有公开实体，包括 double、byte、max_orbitals、
max_radial_points 及全部数组和常量。编译器导出表逐项核对通过，模块元数据标识
确认种类为 MODULE 后排除，不删减任何数据实体。

仅修改过程内 USE，不改变模块头部重导出、变量声明、初始化或表达式。
逐行验证 ONLY 编辑范围，物理行数不变，没有新增超长行。
RMCDHF 文件追加 [state4 快照](../reference/rmcdhf_use_only_state4_snapshot.json)，
绑定 debug manifest；来源链及历史快照篡改负例通过。

Fortitude C121 从 1,950 降至 1,757，默认诊断为 5,813，扩展诊断为 13,245。
差异仅为 193 条 C121，其余诊断按文件、规则、行号和消息比较完全一致。
此前提交的位置清单保留为历史记录，接续工作应重新扫描。

使用本机 Homebrew GNU Fortran 15.3.0、共享 Python 和独立 Release/Check 目录。
两套构建通过；158 个受影响 Release 对象逐字节一致，1,201 份构建模块解压内容一致。
所有修改源码都在当前构建目标中，492 份导入源码的来源链检查通过。
Release/Check 完整 CTest 均 89/89 通过；测试生成探针后全部 1,204 份模块内容也一致。
`git diff --check` 通过。
本地忽略目录 `build-use-only-state4-20261117/` 保存修改前源码、导出列表、哈希、
静态检查、构建及回归日志。本批不新增物理计算，已有科学验证限制保持不变。
