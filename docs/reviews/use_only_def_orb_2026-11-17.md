# USE ONLY 第二十三批：物理参数和轨道状态

日期为 2026-11-17。基准为 0cb432b 加第十二至二十二批工作区修改。

149 份源码补充 181 处 ONLY。非 RMCDHF 源码处理 def_c 117 处，RCI 和跃迁程序处理
orb_c 64 处。按 Release 和 Check 编译命令定位模块，核对公开导出表和源码对符号的
引用。三个跃迁变体继续以 ccmps 引用 speed_of_light_cm_per_s。本批涉及的模块内
导入均位于各自的过程作用域。数值声明、常量和计算语句保持原样。
RMCDHF 及其他程序中余下的 orb_c 导入留待后续批次。

Fortitude 扫描 1,731 份源码，C121 从 779 降至 598。默认诊断共 4,654 条，
扩展诊断共 12,086 条。校正续行产生的行号偏移后，其余诊断的位置、规则和消息一致。

Release 和 Check 构建均通过，CTest 各 89/89 通过。
1,201 个 Release 构建模块接口解压后的内容哈希一致。
在临时副本中对齐行号后，168 个受影响的 Release 编译对象逐字节一致。
RMCDHF 来源校验和 git diff --check 均通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-def-orb-20261117。
