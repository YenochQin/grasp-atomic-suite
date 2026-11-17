# USE ONLY 第十批：PACK 重导出接口

日期：2026-11-17。基准为 bea6655 加未提交的第七至九批修改。本批未提交、未推送。

11 份源码的 `use pack_i` 改为 `use pack_i, only: pack`。
pack_i 将 pack_implementation 的 pack_impl 重导出为 pack，以避开 GNU Fortran 的
内建过程同名绑定问题。本批保留模块重导出实现，只在调用方显式导入既有公开过程；
编译器导出表确认只有 pack，不能误改为直接调用内建 PACK 或更换 ABI。

逐行检查仅追加 ONLY，行数、声明、数值常量和调用不变。
四份 RMCDHF 文件追加 [pack 来源快照](../reference/rmcdhf_use_only_pack_snapshot.json)，
绑定 state4 manifest，历史快照不变；来源链篡改负例通过。

Fortitude C121 从 1,757 降至 1,746，默认诊断 5,802，扩展诊断 13,234。
仅减少 11 条 C121，其余诊断按文件、规则、行号和消息比较完全一致。
此前剩余位置清单是历史记录，后续应重新扫描。

完整 Release 和 Check 构建通过，两套 CTest 均为 89/89 通过。
与第九批基准相比，受影响的 13 个 Release 对象逐字节一致，全部 1,204 个
编译模块接口解压后逐字节一致；PACK 的重导出绑定和接口保持不变。
来源校验通过：492 份导入文件匹配已审核快照，其中 480 份为 Fortran 文件。
`git diff --check` 通过。本批没有修改数值常量、精度声明或计算语句。

使用本机 Homebrew GNU Fortran 15.3.0、共享 Python、独立 Release/Check 目录。
证据在本地忽略目录 build-use-only-pack-20261117，包括修改前源码、导出表、对象及模块
对照、静态诊断和回归日志。未新增物理计算，既有科学验证限制保持不变。
