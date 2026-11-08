# 移除参数不匹配放行选项与接口验证

日期：2026-11-08。基线为 d1260ba 加上一轮未提交的 F2018 标准兼容性修复。
本轮按用户要求移除 `-fallow-argument-mismatch`；过程中先行[移除旧 orbopt 系列](../implemented/orbopt_retirement.md)。

**当前保留目标在 GNU Fortran 15.3.0、`-std=f2018` 下全量编译/链接成功，不含参数不匹配放行选项。**
`-fno-automatic` 仍保留。本轮不是完整的存储状态、全部手写过程接口或物理精度验收。

## 实施内容

1. **MPI 官方显式接口。** `mpi_C` 和 `mpiu.f90` 从 `mpif.h` 改为 `USE mpi`，保留现有整数 handle ABI。
   不增加自定义 `NO_ARG_CHECK` 或重新启用放行选项。MPI 的 choice buffer 本来就允许不同数据类型；
   官方模块提供对应接口，而 datatype/count 等普通参数仍受检查。
   这是 MPI 标准支持的迁移路径；未来 `mpi_f08` 的强类型 handle 改造另行实施。
   [MPI Fortran 绑定说明](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node468.htm)
2. **修正实际句柄和输出约定。** MPI datatype handle 使用默认 INTEGER，不再存入 INTEGER(BYTE)。
   BYTE 仍用于原 CSF 数据元素。为 `startmpi/mpix_bytes` 提供显式接口，并修正初始化例程中
   被错误标成输入的 rank/size/host/length 输出参数。
3. **统一 LP64 BLAS/LAPACK 接口。** 新增
   [suite_linear_algebra](../../src/lib/libmod/linear_algebra_interfaces.f90)，集中声明 DCOPY/DDOT/DSCAL/DSWAP、
   DAXPY/DNRM2/IDAMAX、DGEMV/DGEMM 和 DSPEVX；原 `dspevx_I` 转导同一声明。
   按当前库约定检查默认 INTEGER 为 32 位、DOUBLE 为 64 位。
4. **真实数组区段。** 将指针数组元素作为 BLAS 数组起点的非法用法改为有界的一维区段，
   保持原始起点、元素数与步长；没有通过临时单元素数组丢弃输出，也没有重新解释类型来绕过检查。
   未改合法的 assumed-size/explicit-shape 数组元素序列关联。
5. **精度和读写契约。** LAPACK 的双精度 VL/VU 参数改用双精度常数；INIESTSD 的 BASIS
   实际被初始化/写入，将实现及接口中的错误 `INTENT(IN)` 同步改为 `INTENT(INOUT)`。
   DSPEVX 的 AP 明确为输入输出，其他输出/工作数组按数值库契约声明。
6. **DDOT 返回类型。** 移除错误的整数返回声明，调用方使用统一双精度函数接口。
   在当前环境中，以同一 OpenBLAS 测试 `[1.25,2,3] · [-2,4,1.25]`，数学结果为 9.25；
   旧 INTEGER 声明得到 2，新接口得到正确结果。该示例验证的是 ABI 错误，不是原子计算误差的估计。
   [Netlib DDOT 契约](https://netlib.org/lapack/explore-html/d1/dcc/group__dot_ga2a42ecc597403b22ad786715c739196b.html)

DDOT 修正可能改变首次矩阵求解之后、多根且启用 CI 系数阻尼的正交化分支结果。
这是独立的缺陷修复，不能强求复现旧的错误返回值；本轮没有完成真实多根阻尼计算的物理验证。
默认单根、无该阻尼分支的一轮回归结果见下文。

## 验证结果

环境与上一轮相同：GNU Fortran 15.3.0、Open MPI 5.0.11、系统 OpenBLAS；
`FC` 和 `OMPI_FC` 均为 `/usr/bin/gfortran-15`，MPI/串行/CSFG/OpenMP 开启，Debug 单任务构建。

| 检查 | 结果 |
| --- | --- |
| 当前全量构建 | 成功；3103 个编译条目均使用 GNU 15/F2018，不含 `-fallow-argument-mismatch` 或 `-std=legacy` |
| orbopt 清理 | 四个源码目录、目标与专用驱动已删除；当前目标清单和编译清单无旧目标/目录，活动 bin 中无同名旧程序 |
| 网格补丁器 | 35 个用例通过；suite 当前 21 个网格入口，来源 GRASP 布局未改 |
| BLAS/LAPACK 数值契约 | 9.25 点积、指针区段/负步长/非连续视图拷贝，以及解析 2×2 对称矩阵的本征值、残差和正交性通过 |
| 编译期负向检查 | 合法数组/句柄通过；非法指针元素、错误计数 kind、错误 DDOT 返回声明和 BYTE 句柄均被拒绝 |
| MPI 1/2/4 ranks | BYTE 二维数据、64 位整数、双精度、LOGICAL、CHARACTER 广播及全归约通过 |
| 原版两个串行 solver 一轮 SCF | 共 14 个波函数、混合系数和 MCP 文件与上一轮保留的 gfortran 15 构建逐字节一致 |
| 当前完整 CTest | **36/47 通过，11 项失败，无跳过** |

测试总数由上一轮 44 减去两个随 orbopt 退役的专用测试，再增加五个接口测试，得到 47。
剩余失败名集合仍为原有 NP 九项和 strict restart 两项；没有新增失败。
strict 测试的 baseline 已由旧 orbopt memory 改为原版 `rmcdhf_mem_mpi`，这一设置变化单独记录，
不将相同测试名称误解为整个旧测试环境完全相同。其本轮失败仍发生在期待 `STRICT_CONVERGED` 的断言。

一轮 SCF 对照使用相同 synthetic Z=10、`1s^2`、N=590、H=0.05、ACCY=1.5625e-8 输入，
固定单次外层迭代和单线程 BLAS。预期是等价修复路径输出一致，实际满足；不以该结果宣称 SCF 或活性空间收敛。

## 仍需处理

- NP 摘要空路径和 strict 既有未收敛测试尚未修复。
- 老 Davidson 手写 `_I` 接口仍有 INTENT 属性不一致告警；本轮没有用警告关闭选项隐藏它们。
  全部内部接口统一、`mpi_f08`、大计数/乘法溢出和移除静态存储策略继续按 plan 分项处理。
- LP64 接口不支持把全局整数默认宽度改成 64 位；超出 32 位的 BLAS 计数需另行设计。
- 真实多根、完整 RCI/跃迁参数、跨优化级别和跨编译器物理验证未完成。

本地证据在忽略提交的 `build-gf15-interfaces/audit/`：失败及分批构建日志、
`contract-tests.log`、`ctest-full.log`、`ctest-summary.json`、`first-sweep-parity/comparison.json`、
`legacy-ddot.log`、目标退役清单和源码修改前 hash。验证期间未提交集群作业；构建产物与日志不纳入版本控制。
