# F2018 兼容性与共享接口

更新：2026-11-20，移除旧 kind 接口。文档整理核对基准为 `80f6e76`。本页汇总已实施约束；
逐批审查、编译探查和修改报告由 Git 历史保存。当前构建以
[构建策略](build_profiles.md)及 [Rust 原生输入构建说明](native_csf_descriptors.md)为准。
未完成的数值、状态与平台验收见 [F2018 后续计划](../plans/f2018_migration_gfortran15_2026-11-08.md)。

## 标准、存储与平台接口

- 严格验收使用 GNU Fortran 16，默认告警基线覆盖 16.2 系列。默认生产配置为 `f2018-release`。
  `f2018-check` 和 `f2018-audit` 分别用于运行时及优化后诊断。
  所有配置已移除 legacy、隐式 SAVE 和参数不匹配放行选项。
- 显式数值 kind 直接从 `iso_fortran_env` 导入，例如 `real64`、`int64` 和 `int8`。
  套件已移除 `vast_kind_param` 和其他共享模块重导出的旧 kind 名。
  外部 Fortran 调用方需改用标准名，再重新编译。`grasp_parameters` 导出六个语义容量常量。
- 默认整数保持原有宽度。先以默认整数相乘再转换为 `int64` 仍可能溢出，大计数需单独核查。
- [posix_services](../../src/lib/core/posix_services.f90) 使用 C 互操作连接 POSIX 服务，
  保留目录操作、访问检查及退出状态契约；这是已验证平台范围内的实现，不代表所有系统兼容。
- [array_storage_size](../../src/lib/grasp_shared/array_storage_size.f90) 通过 assumed-rank、
  `STORAGE_SIZE` 和 `SIZE` 计算有效对象的 payload 字节数。
- CSFG、orbopt、strict/diag 变体已退役；历史修复或验收记录不授权恢复这些目标。
  现行范围见[退役说明](rejected_rmcdhf_variants.md)。

## 共享库过程与导入约束

`3e02e0f` 将共享库原 188 份独立接口中的 200 个具体过程迁入实现模块：
`9290` 74 个、`dvd90` 11 个、当时的 `mcp90` 2 个、`rang90` 89 个、`mpiu90` 24 个。
`mcp90` 后续更名为 [`breit_angular_support`](breit_angular_support_migration.md)。

- 实现默认 private，显式公开原 `*_i` 模块 API；已删除的共享库 `*_I.f90` 不再维护副本。
  使用关键字的调用方按当前过程声明传参。原接口中 `itjj::kg1`、`cpath` 路径参数的 INTENT
  已按实际写入行为改为 OUT。
- [count_radial_nodes](../../src/lib/core/radial_node_count.f90) 和
  [pack_csf_subshell_byte](../../src/lib/core/csf_byte_packing.f90) 是各自模块的具体过程。
  接口回归覆盖过程指针、过程实参和关键字调用。
- [mpi_runtime](../../src/lib/mpi/mpi_runtime.f90) 提供无套件状态依赖的
  `initialize_mpi_runtime`、`create_contiguous_mpi_byte_datatype`、`finalize_mpi_runtime`；`mpi_utilities` 保留重导出，避免模块依赖环。
- `davidson_eigensolver::davidson_matrix_vector_operator` 和 `calculate_two_electron_angular_coefficients` 的回调接口仍是契约的一部分；
  测试替身使用相同模块 API，不能通过另一个外部符号替换模块过程。
- Davidson 求解库的现行模块和关键字形参见[命名对照](davidson_eigensolver_migration.md)。
- 这次模块化改变 Fortran 链接符号。依赖旧静态库的外部项目需要重新编译和链接，
  不能将保留源码 API 理解为保留旧二进制 ABI。

`0cb432b`、`5806884` 已完成分批显式导入。
[mpi_state](../../src/lib/mpi/mpi_state.f90) 的模块级 `use mpi` 保留为已复核例外：
它公开转发所用 MPI 实现的符号，不能为清零 C121 而随意收窄接口。
后续调整 ONLY 列表仍须检查重命名、重导出、泛型、回调和条件编译分支。

现行过程名见[共享库过程对照](shared_library_procedure_names.md)。

两个 Biotransform 实现通过 `invert_banded_matrix` 接收 `matrix`、`inverse_workspace`、
`matrix_order`、`determinant`、`singularity_tolerance`、`inversion_status` 和 `leading_dimension`。
矩阵与工作区的显式形状均为 `(leading_dimension, leading_dimension)`。

## 已实施的输入、数值及接口修复

- `80992a9` 修复 GETALDWT 的读取长度、RCI 积分统计、Davidson 接口/状态和 INTERP 边界。
  GETALDWT 使用接口规定的 NCF 读取权重，不能把未定义的旧 NCMIN 随意置零来消除诊断。
  DGEMVMPI 的非单位步长分支在入口已被拒绝；删除它属于不可达代码清理，
  不代表修复了有效输入上的数值错误或新增通用步长支持。
- `9d6ee31` 增加 MCP 读取状态即时检查、MPI 故障终止及 CPATH 分配保护。
- `bea6655` 的 S1/S2 使 RIS4 缺失/不完整缓存安全降级，读取失败按串行/MPI 协议报告；
  JJ2LSJ 唯一标签读取拒绝损坏记录，保留正常输出和选择算法。
- 同批 S3/S4 将限定的 20 个非精确实常量改为 real64，修正两份 BNDINV 接口边界，
  为 JJ2LSJ 指针成员设默认空关联，并报告 MPI worker 输出重定向失败。
  `SLAMCH` 的单精度返回 ABI 保留，不能批量提升所有常量或接口。
- `8cf6b1d`、`80f6e76` 已替换生产 CSF 文本读取入口；
  当前文件、错误传播及低内存契约见[原生描述符说明](native_csf_descriptors.md)。

“未初始化”诊断需要结合赋值路径、有效输入和可达性判断；
诊断条数不等于变量数或已复现缺陷数，不能以统一补零或屏蔽告警代替修复。

## 验证入口与边界

[test/f2018/README.md](../../test/f2018/README.md) 说明强制清单、矩阵、
警告采集和数值回归。配置好独立构建目录后，按改动范围选择已有测试：

```fish
ctest --test-dir build --output-on-failure
../graspkit-tools/.venv/bin/python scripts/check_explicit_contexts.py build \
    --output build/explicit-contexts.json
```

2026-11-17 文档整理未运行新的构建、静态工具扫描或科学计算。
此前模块化的串行/MPI/OpenMP 回归、S1–S4 的 Release/Check 回归及显式导入验证
只证明各自提交和夹具范围；当前原生输入的验收见[独立记录](../reviews/native_csf_descriptor_validation_2026-11-17.md)。

S3 数值对照中，普通 RMCDHF 小夹具两边均执行 6 轮，科学输出一致；
Ni/Zn 初始轨道有浮点末位变化，不能据此宣称所有输出逐字节不变或物理精度提高。
当时 memory 完整 SCF 对照因 macOS 缺少 `/proc/meminfo` 未进入迭代，
该对照需在支持环境补做；也未新增完整 RIS4 物理结果或 MPI SCF 收敛对照。
原实验输入、输出和日志在当时忽略目录 `build-static-numeric-20261117/`，
其中 `numerical-comparison.json` 汇总结果；这些本地产物不是仓库可移植依赖。

历史告警预算覆盖缺口仍单列，不能以 CTest 通过代替告警门禁或物理验收。
后续事项见[静态诊断清单](../plans/static_diagnostics_remaining_2026-11-17.md)。
[警告基线](../reference/f2018/gfortran16-warning-baseline.md)登记风险，不将已有诊断视为安全豁免。

## 来源校验

[原始 RMCDHF 清单](../reference/rmcdhf_upstream_snapshot.json)及
[兼容性清单](../reference/rmcdhf_f2018_compatibility.json)和后续快照链继续由
[check_upstream.py](../../test/npol/check_upstream.py)使用。
删除说明性报告不删除这些机器可读来源记录，也不放宽精确范围或哈希检查。
190 份受此次 kind 清理影响的 RMCDHF 源文件登记在 [kind 清理快照](../reference/rmcdhf_kind_alias_cleanup_snapshot.json)。
新增源码修改须按现行校验链登记，不能覆盖旧 manifest。
