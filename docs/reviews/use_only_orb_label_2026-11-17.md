# USE ONLY 第二十五批：轨道状态和角动量标签

日期为 2026-11-17。基准为 0cb432b 加第十二至二十四批工作区修改。

174 份源码处理 180 处缺少 ONLY 的导入，包括 150 处 orb_c、28 处 jlabl_c
和 2 处 biorb_c。范围覆盖四个普通 RMCDHF 变体剩余的 64 处 orb_c。
对照 Release 和 Check 的模块导出表及源码引用，补充 176 处 ONLY，保留
ncftot、ncfblock、labj 和 labp 的原有映射。两份 CMake 未使用的 orbord.f90
各移除两条未引用的导入。复查 src 和 test，三模块的其余导入均含 ONLY。
72 份 RMCDHF 改动加入 [来源快照](../reference/rmcdhf_use_only_orb_label_snapshot.json)。
来源校验及前驱快照篡改负例均通过。

Fortitude 扫描 1,731 份源码，C121 从 418 降至 238。默认诊断共 4,294 条，
扩展诊断共 11,726 条。校正续行及删行造成的行号偏移后，其余诊断的位置、规则和消息一致。

Release 和 Check 构建均通过，CTest 各 89/89 通过。
1,201 个 Release 模块接口解压后的内容哈希一致。在临时副本中对齐源码行号后，
207 个主构建的 Release 编译对象逐字节一致。两份 orbord.f90 使用相应串行和
MPI 程序的编译参数单独对照，Release 对象逐字节一致，Check 对象去除调试信息后
逐字节一致。git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-orb-label-20261117。
