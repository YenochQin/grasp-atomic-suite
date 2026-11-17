# USE ONLY 第十九批：径向网格和波函数模块

日期为 2026-11-17。基准为 0cb432b 加第十二至十八批工作区修改。

106 份源码补充 152 处 ONLY，其中 grid_c 87 处，wave_c 65 处。
按基准 Release 编译器导出表保留全部公开符号，包括 double 精度别名和
radial_workspace_capacity、max_orbitals 容量参数。源码对照确认改动限于导入语句。
本批范围为已纳入 CMake 构建的非 RMCDHF 外部程序单元。RMCDHF 和模块内剩余导入待处理。

Fortitude 扫描 1,731 份源码。C121 从 1,429 降至 1,277，默认诊断共 5,333 条，
扩展诊断共 12,765 条。校正续行造成的行号偏移后，其余诊断的位置、规则和消息一致。

复用已有的 Release 和 Check 构建目录，两套构建均通过，CTest 各 89/89 通过。
1,201 个 Release 生产模块接口的解压内容哈希一致。在临时副本中对齐源码行号后，
122 个受影响编译对象逐字节一致。RMCDHF 来源校验和 git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-radial-20261117。
