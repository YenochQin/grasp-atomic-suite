# USE ONLY 第十六批：跃迁振子强度状态模块

日期为 2026-11-17。基准为 0cb432b 加第十二至十五批工作区修改。

在 rtransition、rtransition_mpi 和 rtransition_phase 的 21 份源码中为
osc_c 补充 21 处 ONLY，保留全部 25 个公开符号，包括 double 精度别名。
对照模块导出表和折叠续行后的源码，确认改动限于导入语句。

Fortitude 默认诊断 5,664 条，扩展诊断 13,096 条。C121 从 1,629 降至 1,608。
按文件、规则和消息比较，其他诊断一致。续行使后续行号顺延。

清除 80 个历史 CMake 构建树前，记录上一批 1,204 个 Release 模块文件的
解压内容哈希。重新构建和回归后，模块文件清单及哈希全部一致。
在临时副本中对齐修改前后的源码行号，使用相同 Release 编译参数得到的
21 个对象逐字节一致。来源校验和 git diff --check 均通过。

Release 和 Check 全新构建均通过，CTest 各 89/89 通过。使用本机 Homebrew
GNU Fortran 15 和共享 Python。当前验证证据位于忽略目录
build-use-only-osc-20261117。历史构建树内的独立审计及计算记录移至
build-history-evidence，同级日志和静态扫描记录保留。
