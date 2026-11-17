# USE ONLY 第十七批：跃迁程序的 JJ2LSJ 数据导入

日期为 2026-11-17。基准为 0cb432b 加第十二至十六批工作区修改。

rtransition、rtransition_mpi 和 rtransition_phase 各有 5 份源码为
jj2lsjbio_c 补充 ONLY。按基准 Release 编译器导出表保留全部 15 个公开符号，
包括 double 精度别名。源码对照确认改动限于这 15 条导入语句及其续行。

Fortitude 扫描 1,731 份源码。C121 从 1,608 降至 1,593，默认诊断共 5,649 条，
扩展诊断共 13,081 条。校正续行引起的行号偏移后，其余诊断的位置、规则和消息一致。

从空目录构建基准 Release，修改后重建 Release，并从空目录构建 Check。
两套 CTest 各 89/89 通过。生产目录内 1,201 个 Release 模块接口的解压内容哈希一致。
在临时副本中对齐源码行号后，15 个受影响编译对象逐字节一致。
RMCDHF 来源校验和 git diff --check 均通过。

本机使用 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
验证日志和对照数据位于忽略目录 build-use-only-jj2lsjbio-20261117。
