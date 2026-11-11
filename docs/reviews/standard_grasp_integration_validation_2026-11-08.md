# 普通 GRASP 导入验证记录

> 2026-11-11 状态更新：文中 strict/diag 独立求解器已按用户决定[移除](../implemented/rejected_rmcdhf_variants.md)；相关目标、命令和测试数量仅作对应版本的历史记录。

日期：2026-11-08。suite 修改前基线：`179ca9ecfd96d54a5617d811b888299a611c13b3`。
变更范围见[实施说明](../implemented/standard_grasp_integration.md)。
本记录区分本次接入验证、既有失败与未执行的物理验收。

## 验证环境与结果

GNU Fortran 15.3.0、CMake 4.4.5、系统 OpenBLAS、Open MPI 5.0.11；Debug 构建使用原有
`-fno-automatic -fallow-argument-mismatch` 兼容选项，不是严格 F2018 模式。
配置时指定 `gfortran-15`，测试中的编译探针也使用该版本。MPI 的系统模块/库在本机
可以完成构建及下述通过的运行测试，不据此推广为所有 GCC 15/MPI 组合均兼容。

| 项目 | 结果 |
| --- | --- |
| MPI/串行/CSFG/OpenMP 全部启用的 suite 构建 | 完成编译和链接；包含本次十个目标 |
| MPI、CSFG、OpenMP 关闭，串行应用开启 | 六个新串行目标构建通过，未创建新 MPI 目标 |
| 网格补丁器 | 35 个用例通过，含真实新增入口的修改/逐字节恢复 |
| 默认容量检查 | `--layout atomic-suite --nnnp 590 --check` 无待修改文件 |
| 非默认容量 | 共享数组和普通波函数接口在 790/800/149 配置下通过 |
| 非默认容量反例 | 临时恢复原 `tail_I.f90` 后，测试按预期在过程指针关联处因数组形状不一致而编译失败 |
| 全部 CTest，合并受限环境首轮与本机 MPI 重跑结果 | **28/39 通过，11/39 失败**；不是全绿验收 |
| 11 个失败的原基线对照 | 未修改 suite 基线及原版准备工具下 **11/11 同样失败** |

首轮全量 `--parallel 4` 构建在既有 `gj90/gj90_mpi` 同时编译 `gethfd_I.f90` 时出现
`Cannot rename module file gethfd_i.mod0`。这两个目标原本共享模块输出目录，相关 CMake
未在本次改动。随后使用 `--parallel 1` 完成全量构建；新串行目标的独立配置可用
`--parallel 4` 构建。该既有并行模块输出冲突没有在本次五类程序接入中修复。

## CTest 的已通过范围

已通过项目包含共享容量、QUAD、严格停止逻辑、轨道事务/轮次逻辑、P0/START 单元检查、
新准备工具参与的 `p0_trace_mpi_1/2`、密度串行及 1/2/4-rank 回归、jj2lsj 线程/稳定排序、
普通和 CSFG 容量/文件格式，以及 NP 的 kernel/bounds/input/来源保护/预期 MPI 错误路径。

沙箱首轮禁止 MPI 的本机 socket，出现 `pmix_ifinit: socket() failed with errno=1`。
这些环境阻断与实际计算失败分开记录；在允许本机 MPI 通信的环境重跑全部首轮失败项后，
P0 链、密度 MPI 和 NP 预期错误路径通过，留下下列 11 项。

## 已有失败的独立复现

从上述 Git 提交归档未修改 suite，在独立目录用同一 GCC 15 Debug 设置重新构建。
另外从工作区原始 `grasp/` 制作隔离副本，只将该副本 NNNW=100 对齐为 suite 的 127，
用 GCC 15 构建原版 `rwfnestimate/rangular_mpi`。原工作区来源没有修改。
原基线使用这两个原版准备工具运行相同失败测试；未使用本次导入后的应用实现。

| 失败组 | 测试 | 原基线与当前结果 |
| --- | --- | --- |
| NP 摘要输出，9 项 | `npol_code_compose`、`npol_code_output`、`npol_code_entrypoints`、四个 `npol_code_consumers_rmcdhf*`、`npol_code_mpi_broadcast_1/2` | 都报告摘要路径为空，`cannot open summary: Cannot open file ''`；MPI 正常广播路径也因此失败 |
| strict 重启，2 项 | `strict_mem_restart_mpi_1/2` | 固定样例重启后没有得到测试要求的 `SCF termination: STRICT_CONVERGED`，断言失败 |

这证明在本机工具链及这些夹具下，失败不依赖本次导入和参数归属修改。
没有把 NP 的问题直接归因于编译器缺陷，也没有完成 strict 不收敛的根因分析。
本次未修改 NP、SCF、停止判据或验收阈值来使测试变绿；保留这些失败以供后续单独处理。

## 来源和数据边界

815 个导入文件中，只有 `rwfnestimate90` 的六个 Fortran 文件因容量导入而变化；
其他导入文件与[来源清单](../reference/standard_grasp_upstream_snapshot.json)逐字节匹配。
来源 `grasp/` 的对应文件全部保持原 hash。新程序的 CMake 源文件清单未擅自扩充。

执行的是工程检查和已有小型回归夹具，没有提交 SLURM 作业，没有开展真实原子的
活性空间收敛或性能研究。RCI、双正交和跃迁目标已验证构建，但本次没有完成它们的
新旧完整数值对照；不能据此声明新计算链的物理准确性。

本地日志和机器可读汇总保存在忽略提交的
`build-gf15-standard-import/validation/`，包括配置、构建、首轮 CTest、本机 MPI 重跑、
未修改基线复现、来源适配检查及容量反例。该目录不是仓库发布内容。
