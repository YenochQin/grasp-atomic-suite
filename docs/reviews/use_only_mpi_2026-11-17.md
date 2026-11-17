# USE ONLY 第二十六批：MPI 和串行模块导入

日期为 2026-11-17。基准为 0cb432b 加第十二至二十五批工作区修改。

143 份文件处理 166 处缺少 ONLY 的导入。覆盖 mpi_c 95 处、mpi_s 27 处、
直接导入 mpi 16 处，另有 28 处来自 rmcdhf_mem_c、wchblk_c、tatb_c、
vpilst_c、hmat_c 和 vinlst_c。保留所有已有重命名。mpix_startup 不直接
引用 MPI 符号，移除一条未使用的 USE。

逐一核对 Release 和 Check 的模块搜索顺序、编译器导出表以及测试源文件的
条件编译分支。mpi_c 中模块级的 use mpi 仍负责重导出 MPI 的完整公开接口，
保留该语句以维持现有接口和对不同 MPI 实现的兼容性。
74 份 RMCDHF 改动加入 [来源快照](../reference/rmcdhf_use_only_mpi_snapshot.json)。
来源校验和前驱快照篡改负例均通过。

Fortitude 扫描 1,731 份源码，C121 从 238 降至 72。默认诊断共 4,128 条，
扩展诊断共 11,560 条。校正续行和删行的行号偏移后，其余诊断的位置、规则
和消息一致。MPI 相关的剩余 C121 为 mpi_c 中上述公开重导出。

Release 和 Check 构建均通过，CTest 各 89/89 通过。
1,201 个 Release 构建模块接口解压后的内容哈希一致。
在临时副本中对齐源码行号后，169 个受影响的 Release 编译对象逐字节一致。
git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 5.0.11 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-mpi-20261117。
