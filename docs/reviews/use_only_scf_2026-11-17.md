# USE ONLY 第二十一批：SCF、哈密顿矩阵和过程接口

日期为 2026-11-17。基准为 0cb432b 加第十二至二十批工作区修改。

116 份源码处理 153 处缺少 ONLY 的导入，覆盖 scf_c 80 处、hmat_c 31 处、
itjpo_i、ispar_i、iq_i、jqs_i、dspevx_i、gdvd_i、count_i 和 nsmdat_c 共 42 处。
其中 151 处明确列出原有公开符号。JJ2LSJ 的 setls_equivalent_csfs 删除 2 处重复导入，
通过宿主关联继续使用外层 setls 从同一模块导入的 itjpo 和 iq。
按各构建目标的 Release 模块导出表核对，保留精度别名、容量参数和过程接口。
memory_man 的生产模块与测试替身公开不同符号，留待按目标分开核对。

四个普通 RMCDHF 变体的 96 份改动加入
[来源快照](../reference/rmcdhf_use_only_scf_snapshot.json)，衔接径向数据快照。
来源校验与前驱快照篡改负例均通过。

Fortitude 对 1,731 份源码的 C121 从 1,114 降至 961。默认诊断共 5,017 条，
扩展诊断共 12,449 条。校正增删行引起的位置偏移后，其余诊断一致。

Release 和 Check 构建均通过，CTest 各 89/89 通过。1,201 个 Release 生产模块接口
解压后的内容哈希一致。在临时副本中对齐源码行号后，129 个受影响的编译对象
逐字节一致。git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-scf-20261117。
