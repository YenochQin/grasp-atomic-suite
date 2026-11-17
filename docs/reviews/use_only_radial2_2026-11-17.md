# USE ONLY 第二十批：RMCDHF 和模块内的径向数据导入

日期为 2026-11-17。基准为 0cb432b 加第十二至十九批工作区修改。

107 份文件补充 163 处 ONLY，其中 grid_c 87 处、wave_c 76 处。范围包括
四个普通 RMCDHF 变体的 88 份源码、18 份模块源码和 lib9290_quad 测试。
参照基准 Release 导出表保留全部公开符号，包括精度别名、容量参数及测试中的
八个既有重命名。在 src 和 test 中复查，两模块的导入均已包含 ONLY。
RMCDHF 的 88 份改动加入 [来源快照](../reference/rmcdhf_use_only_radial_snapshot.json)，
衔接 EIGV 快照。来源校验与前驱快照篡改负例均通过。

Fortitude 对 1,731 份源码的 C121 从 1,277 降至 1,114。默认诊断共 5,170 条，
扩展诊断共 12,602 条。校正续行造成的行号偏移后，其余诊断的位置、规则和消息一致。

Release 和 Check 构建均通过，CTest 各 89/89 通过。1,201 个 Release 生产模块接口
解压后的内容哈希一致。用临时副本对齐行号后，116 个生产编译对象逐字节一致。
积分测试的 Release 对象逐字节一致，Check 对象去除调试信息后逐字节一致。
git diff --check 通过。

使用本机 Homebrew GNU Fortran 15.3.0、Open MPI 和共享 Python。
本批日志及对照数据位于忽略目录 build-use-only-radial2-20261117。
