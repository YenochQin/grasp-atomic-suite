# USE ONLY 第二十二批：分配模块与测试替身

日期为 2026-11-17。基准为 0cb432b 加第十二至二十一批工作区修改。

182 份源码为 memory_man 补充 182 处 ONLY，覆盖 src 和 test 中此前缺少
ONLY 的全部 memory_man 导入。根据 Release 和 Check 的编译命令逐目标核对
模块搜索顺序、编译器导出表和源码引用。每个列表使用对应构建目标共同提供的符号。

串行和 MPI 的 setmcp_mem 各由生产、常规输入测试和故障注入测试目标编译。
这两份共享源码均只导入泛型 alloc。检查两套构建的 12 个相关对象符号，确认
生产和常规输入测试引用生产分配过程，故障注入测试引用各自的测试替身过程。
44 份 RMCDHF 源码改动加入 [来源快照](../reference/rmcdhf_use_only_memory_snapshot.json)。
来源校验和前驱快照篡改负例均通过。

Fortitude 对 1,731 份源码的 C121 从 961 降至 779。默认诊断共 4,835 条，
扩展诊断共 12,267 条。其余诊断的文件、行号、规则和消息一致。

Release 和 Check 构建均通过，CTest 各 89/89 通过，其中包括串行和 MPI 的
分配故障注入测试。1,201 个 Release 模块文件的解压内容哈希一致。
按原目标编译参数对照，206 个受影响的 Release 对象逐字节一致，包含两种
测试替身的调用对象。git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-memory-20261117。
