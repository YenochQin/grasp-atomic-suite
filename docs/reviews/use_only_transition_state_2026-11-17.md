# USE ONLY 第十八批：轨道变换、分块和本征解模块

日期为 2026-11-17。基准为 0cb432b 加第十二至十七批工作区修改。

138 份源码补充 164 处 ONLY，其中 biorb_c 50 处、blk_c 21 处、eigv_c 93 处。
按基准 Release 编译器导出表保留全部公开符号，包括 double 和 max_orbitals。
对照确认源码改动限于导入语句。串行和 MPI 的两份 orbord.f90 未列入当前 CMake
构建，保留待处理。四个 RMCDHF 变体共 28 份改动加入
[来源快照](../reference/rmcdhf_use_only_eigv_snapshot.json)，衔接 POS 快照。

Fortitude 扫描 1,731 份源码。C121 从 1,593 降至 1,429，默认诊断共 5,485 条，
扩展诊断共 12,917 条。校正续行引起的行号偏移后，其余诊断的位置、规则和消息一致。

复用上一批的 Release 和 Check 构建目录，两套构建及 CTest 均通过，测试各 89/89。
1,201 个 Release 生产模块接口的解压内容哈希一致。在临时副本里对齐行号后，
154 个受影响编译对象逐字节一致。RMCDHF 来源校验、快照篡改负例和
git diff --check 均通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-transition-state-20261117。
