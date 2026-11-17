# USE ONLY 第二十四批：RMCDHF 物理参数和轨道状态

日期为 2026-11-17。基准为 0cb432b 加第十二至二十三批工作区修改。

四个普通 RMCDHF 变体的 136 份源码补充 180 处 ONLY。处理全部 116 处
剩余 def_c 导入，并对四个变体各选相同的 16 个例程，处理 64 处 orb_c 导入。
根据 Release 和 Check 编译命令核对模块公开导出表和源码引用。保留
rwtdum => wt、iiqa => csf_subshell_occupations、ncftot => csf_count 和
ncfblock => csf_count。四个变体仍有 64 处 orb_c 导入待处理。
136 份改动加入 [来源快照](../reference/rmcdhf_use_only_def_orb_snapshot.json)。
来源校验和前驱快照篡改负例均通过。

Fortitude 扫描 1,731 份源码，C121 从 598 降至 418。默认诊断共 4,474 条，
扩展诊断共 11,906 条。校正续行的行号偏移后，其余诊断的位置、规则和消息一致。

Release 和 Check 构建均通过，CTest 各 89/89 通过。
1,201 个 Release 模块接口解压后的内容哈希一致。在临时副本中对齐源码行号后，
147 个受影响的 Release 编译对象逐字节一致。git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-rmcdhf-def-orb-20261117。
