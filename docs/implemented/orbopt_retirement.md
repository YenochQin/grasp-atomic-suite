# 移除旧 rmcdhf_orbopt 程序

日期：2026-11-08。依据用户决定，删除旧轨道优化系列的四个实现及构建/安装入口：

| 已移除目标 | 已移除源码目录 |
| --- | --- |
| `rmcdhf_orbopt` | `src/appl/rmcdhf90` |
| `rmcdhf_orbopt_mem` | `src/appl/rmcdhf90_mem` |
| `rmcdhf_orbopt_mpi` | `src/appl/rmcdhf90_mpi` |
| `rmcdhf_orbopt_mem_mpi` | `src/appl/rmcdhf90_mem_mpi` |

当前保留原版 `rmcdhf`、`rmcdhf_mem`、`rmcdhf_mpi`、`rmcdhf_mem_mpi`，以及独立
`rmcdhf_mem_strict_mpi`、`rmcdhf_mem_diag_mpi` 和 `rmcdhf_csfg_mpi`。
共享数值库、统一参数文件和保留求解器的计算流程不因退役操作改变。

## 配套清理

- 删除旧 orbopt 的专用策略测试、运行器、SLURM 模板及测试夹具；不删除外部计算数据。
- 共有的波函数/能量比较和 strict/diag 测试迁到 [test/rmcdhf](../../test/rmcdhf/README.md)。
  strict 重启测试的对照目标改为原版 `rmcdhf_mem_mpi`，未放宽其收敛断言。
- 删除 `rmcdhf_pair_transaction_logic`、`rmcdhf_orbopt_round_state_logic` 两个已无对应实现的测试。
- F2018 固定一轮差分脚本只覆盖保留的两个串行原版 solver；此前包含 orbopt 的四版本记录保留为历史证据。
- `atomic-suite` 网格布局删去四个旧入口，当前为 **21 个网格入口、4 个 RCI 续算入口**。
  `grasp2018` 布局及来源清单里的上游 `rmcdhf90*` 路径保持原义，不作全局替换。
- 清理 suite 的 `bin/` 和当前接口修复 build 中同名旧可执行文件。
  历史 build 基线保留用于复现，不代表现行可构建目标。

旧 [README](../archive/orbopt/README.md) 和 [实验记录](../archive/orbopt/RESULTS.md) 已归档。
被移除的源码和运行脚本可在 Git 提交 `d1260ba` 追溯；没有将其复制到新的生产目录或提供替代命令别名。

此操作先于移除参数不匹配放行选项的最终回归。后续验证应使用更新后的目标和测试清单，
不能直接以旧的测试总数判断通过率；历史实验结果也不构成当前物理验证。
