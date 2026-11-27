# Suite 原生读取 CSF 计算描述符 V3 计划

2026-11-28 归档。核心接入以 `a2ce6c6` 为源码基线，于 2026-11-21 提交为 `30e35024`。
`612ce512` 补齐高 J 及计算边界，`9d0a85fd` 加入 V3 块索引和文件双向校验。
以 Suite `5e5cfe6d` 复验，40/40 项专项 CTest 通过。本页保留原接入设计和验收用例。
格式主规范由 rCSFs 仓库的 `docs/CSF_DESCRIPTOR_V3_PLAN.md` 维护。
[当前生产输入](../implemented/design-0001-native-csf-descriptors.md) · [专项验收记录](../reviews/evidence-0014-native-csf-descriptor-validation.md#v3-专项验收2026-11-28)

## 输入和内存索引

M 为 header 的有序 peel 子壳层数。每条 CSF 有 `4M+2` 个非空逻辑 Int8，按轨道重复 `(qᵢ, vᵢ, 2Jᵢ+1, 2Kᵢ+1)`。最后两列依次是正整数 `total_two_j_plus_one` 和独立的 `parity`（+1 或 −1）。所有角动量码取 1 至 127。初始累计码为 1；空壳或满壳为 `(qᵢ,0,1,前一步累计码)`，其中 qᵢ 分别取 0 或容量。M>0 时最后的累计码等于总角动量码。M=0 时纯闭核行是 `[1,+1]`。

| 计算需要的量 | Suite 侧实施位置 |
| --- | --- |
| 占据、物理 v、子壳层 2Jᵢ+1 | 按 core、peel 顺序安装轨道表，将 q 填入 `csf_subshell_occupations`，将 `(v,0,2Jᵢ+1)` 填入 `csf_subshell_quantum_codes`。共同闭核满占据、v=0、自身及累计码均为 1。普通 Suite 支持的附加态标签 w 为 0。 |
| 每个轨道结束处的 2Kᵢ+1 | 建立按**轨道索引、CSF 索引**访问的累计数组。现有 `csf_coupling_quantum_codes` 及其访问器第一维按开壳层耦合步编号。为轨道索引提供独立的明确访问入口，逐个迁移使用者。 |
| 总 J 与宇称 | 两列各自读取和验证。总 J 查询返回正的 2J+1，宇称查询返回独立符号。仍要求旧内核有符号末槽的适配点，可在内存中构造 `parity * total_two_j_plus_one`。 |

## 实施顺序

1. **Rust reader 与 C/Fortran 边界。** 在 Suite 自有的 Rust reader 中登记版本 3，核对列语义和次序、header 字节哈希、轨道顺序、行数、显式块边界和可选零阶计数。按宽整数解码物理 2J、2K（码值减一）并验证电子数、占据容量、宇称、子壳层态、每一步的三角条件/步长和终点。重用现有有界批次、顺序读块及错误处理。使 Fortran 可明确查询版本；若扩充 C 接口，同步头文件、ABI 版本策略与绑定。生产路径要求 `descriptor_version=3`，同宽的 V2 文件在打开时拒绝。[Rust 传输契约](../../src/lib/csf_descriptor/src/descriptor_transport.rs) · [Rust reader](../../src/lib/csf_descriptor/src/descriptor_reader.rs) · [C 接口](../../src/lib/csf_descriptor/include/rcsfs_descriptor.h) · [Fortran reader](../../src/lib/core/csf_descriptor_reader.f90)
2. **rangular 与 RCI。** 在块读入时安装本块的占据、子壳层码、逐轨道累计码和总对称性。对每个 bra/ket，继续按占据构造双方共同满壳及其余参与轨道的有序列表。若参与轨道为 p₁…pₛ，则第 t 个中间工作槽直接取各侧在 pₜ₊₁ 处的累计码（t=1…s−1）。s≤1 时没有有效中间槽；单壳层算符例程对总 J 或临时活动壳层的单独处理仍须核对。共同轨道列表、算符选择规则、MCP 角系数和径向积分继续由相应计算例程产生。[成对准备](../../src/lib/core/csf_pair_quantum_numbers.f90) · [量子数查询](../../src/lib/core/csf_quantum_queries.f90) · [RCI 读块](../../src/lib/core/rci_csf_block_reader.f90)
3. **RMCDHF 与 MPI。** 普通、内存及 MPI 版本均保持 MCP 先行，核对闭核数、轨道数、电子数、CSF 总数、显式块长度/标签与进程配置。逐行验证完整输入，长期保留全基组占据和块对称性；沿用 MCP 提供角系数。MPI 由 root 读取和验证，按消费程序广播占据或完整块，并传播错误状态。zero-first 作用于各块有序行的前缀。[RMCDHF 占据路径](../../src/lib/core/rmcdhf_csf_input.f90) · [MPI 块广播](../../src/lib/mpi/mpi_csf_block_input.f90) · [MPI 占据广播](../../src/lib/mpi/mpi_rmcdhf_csf_input.f90)
4. **其他普通 CSF 消费者。** 核查 JJ2LSJ、性质程序、跃迁及双正交变换中直接按步访问、复制或输出耦合码的调用点。全基组 reader 剔除全空轨道时，逐轨道累计数组须与轨道表同步压缩。左右基组合并并插入空轨道时，要分别继承各侧在前一轨道的累计码；继续保留左右轨道的相对次序、两侧行号与块前缀。按步的遗留入口可由开壳层位置映射到对应轨道。header-only 的径向初猜只消费轨道表和电子数。[全基组读入](../../src/lib/core/csf_full_basis_reader.f90) · [双基组合并](../../src/lib/core/csf_basis_merging.f90) · [JJ2LSJ](../../src/apps/jj2lsj/jj2lsj_conversion.f90)

## 验收

- 构造有解析预期的小基组。对 `1s² 2s¹ 3s¹ 4s¹` 在 3s 取 K=0 和 K=1 的两条 CSF，完整累计码分别为 `[1,2,1,2]`、`[1,2,3,2]`。对 B=`1s² 2s² 3s¹ 4s⁰`，完整累计码为 `[1,1,2,2]`。三条都以 `[2,+1]` 结尾。逐项断言 A₀/B 的成对工作区为 `[1,2]`、`[2,2]`；交换 bra/ket 后结果随之交换；B/B 只有一个参与轨道。
- 覆盖纯闭核、奇宇称、部分占据而自身 Jᵢ=0、相同 Jᵢ 不同 v、相邻同对称性块、空轨道插入和剔除。拒绝 q/v/角动量非法组合、末两列颠倒、宇称与轨道不符、累计码未正确继承、总 J 不闭合和损坏的文件对。检查读入失败后的句柄和数组状态。[原生 reader 测试入口](../../test/test_native_descriptors.py) · [双基组测试](../../test/test_descriptor_merge.py)
- 串行及多 rank 检查块边界、局部/全局 CSF 编号、每侧有效耦合码和 zero-first 截断点；用小型已知基组检查选定的角系数、MCP 结构、RCI 矩阵元和 RMCDHF 输入契约。数值验收记录输入、态匹配和容差；仅编译成功不能作为这些结果的证据。
