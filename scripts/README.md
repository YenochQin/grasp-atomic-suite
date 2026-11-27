# GRASP 辅助脚本

2026-11-16：网格工具同时识别原版及语义变量名，命令行/TOML 参数保持不变。
`atomic-suite` 覆盖十五个普通网格入口和两个 RCI 续算入口。
集中配置及 `grasp_parameters` 转导的正式名称见[共同参数说明](../docs/implemented/design-0003-shared-grid-parameters.md)。

大规模 RMCDHF MPI 停顿采样使用 [sample_rmcdhf_mpi.py](sample_rmcdhf_mpi.py)，
只读附着到指定 Linux PID 树，记录线程状态、等待位置、I/O、缺页和系统压力。
用法、原版流程对照及验证边界见 [RMCDHF MPI 诊断](../docs/implemented/design-0021-rmcdhf-mpi-diagnostics.md)。

当前生产 CSF 输入与构建方式见[原生描述符说明](../docs/implemented/design-0001-native-csf-descriptors.md)。
计算程序直接链接 Rust reader，不调用下列开发/对照脚本。

- [compare_rangular_versions.py](compare_rangular_versions.py) 比较默认和显式描述符路径及不同程序版本；
  `--descriptors` 与 `--header` 配套提供。历史参考程序可显式选择 `--reference-text-oracle`；
  两个历史程序的离线文本对照使用 `--legacy-text-oracle`。
- [compare_rangular_mpi_inputs.py](compare_rangular_mpi_inputs.py) 对照同一 MPI 程序的默认/显式原生输入，
  要求 `--descriptors`、`--header`；保留逐 rank 计时和 MCP 哈希。
- [verify_csf_descriptors.py](verify_csf_descriptors.py) 调用不安装的测试 oracle，核对原文本与原生数组。
- [旧 GCSDESC1 转换器](archive/prepare_csf_descriptors.py) 已归档，仅服务历史格式夹具；当前程序拒绝其输出。

角系数文件的解析、统计和原版/候选程序对比见
[MCP 解析与比较工具](../docs/implemented/design-0024-mcp-compare.md)，入口为
[mcp_compare.py](mcp_compare.py)。支持普通串行与 MPI 输出，不运行原子计算。

## 批量修改 grasp-atomic-suite 径向网格

`patch_grasp_grid.py` 根据脚本文件的实际位置找到与 `scripts/` 同级的 `src/`。
要修改某个 Suite 副本，请运行该副本 `scripts/` 中的脚本。命令行 `--nnnp` 和 TOML
`nnnp` 都设置 [cal_constants.f90](../src/lib/grasp_shared/cal_constants.f90)
中的编译期容量 `max_radial_points`。工作区容量按 `max_radial_points + 10` 自动派生。
只传容量时，脚本只修改该中央定义。有限核默认点数引用 `max_radial_points`，因此会随容量改变。
如需指定有限核运行时的默认点数，使用 `--n` 或 TOML `n`。

脚本使用 Python 3.11 或以上版本的标准库。本工作区使用 `graspkit-tools/.venv`。
每组参数使用独立的 Suite 源码副本。修改完成后，按[构建说明](../README_ZH.md)全量重编译
所有实际使用的串行和 MPI 程序。交互输入和旧 RCI `.res` 可覆盖新默认网格，
`GRASP_COUNT_ACCY` 可覆盖工作容差。数值收敛须通过独立计算验证。
见[容量和默认值](../docs/implemented/design-0003-shared-grid-parameters.md)及[网格实验方案](../docs/plans/plan-0001-radial-grid-validation.md)。

### TOML 配置

将 [config.example.toml](config.example.toml) 复制为脚本同目录的 `config.toml`。
备份和恢复的相对路径以配置所在目录为基准，也可使用绝对路径或 `~/`。
本地配置由 Git 忽略。

```toml
RUN_MODE = "preview"
PRINT_DIFF = true

[GRID_PARAMETERS]
nnnp = 1990
```

`GRASP_SOURCE` 和 `SOURCE_LAYOUT` 已移除。使用旧配置时，删除这两个字段。
`RUN_MODE` 默认为 `"preview"`，`PRINT_DIFF` 默认为 `true`。
正常修改需在 `[GRID_PARAMETERS]` 中填写至少一个参数。省略的网格参数保留原赋值。

```fish
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py
```

先预览改动，再设置 `RUN_MODE = "apply"` 写入，最后可用 `RUN_MODE = "check"` 验证。
检查一致时退出 0，需要修改时退出 1，配置错误或源码布局不匹配时退出 2。

### 命令行

显式命令行请求独立于 `config.toml`。以下命令预览当前脚本所属 Suite 的容量修改。

```fish
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py --nnnp 1990 --diff
```

在同一命令后加 `--apply` 写入，或加 `--check` 检查。以下选项控制各入口的默认值。

| CLI / TOML | 默认值的作用范围 |
| --- | --- |
| `--n` / `n`、`--h` / `h`、`--rnt-scale` / `rnt_scale` | 有限核分支和 RCI 的公共初始化。尺度系数除以核电荷后得到 `radial_grid_scale`。 |
| `--point-n` / `point_n`、`--point-h` / `point_h`、`--point-rnt-scale` / `point_rnt_scale` | 有点核分支的入口。 |
| `--hp` / `hp` | 两种核模型的公共渐近步长。设为 0 时选用指数网格。 |
| `--accy` / `accy` | 工作容差，包括两个 RCI 续算入口。设为 0 时使用 `transformed_grid_step**6`。 |

脚本检查十五个网格入口、串行与 MPI 共用的 RCI 续算容差赋值和六项共享容量的转导，并扫描整个 `src/`，
发现未知网格初始化或重复的容量定义时退出。显式设置 H（包括点核 H）或 ACCY 时，
脚本在 RMCDHF 的 ACCY 交互提示前同步容差。

### 扩容并固定有限核默认点数

使用 [config.capacity-only.example.toml](config.capacity-only.example.toml) 中的
`CAPACITY_ONLY = true`，或传入 `--capacity-only`。该模式要求同时提供 `nnnp` 和 `n`，
且只接受这两个网格参数。下面以扩容到 1990 并保留 590 个有限核默认点为例。

```fish
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --capacity-only --nnnp 1990 --n 590 --diff
```

### 备份与恢复

应用补丁时，可通过 `--backup-dir /path/to/new-backup` 或 TOML 顶层
`BACKUP_DIRECTORY` 指定尚不存在的备份目录。脚本保存改动文件的原始字节和校验和。
预览、检查及无变化的重复应用不会创建备份。

```fish
../graspkit-tools/.venv/bin/python scripts/patch_grasp_grid.py \
  --restore /path/to/new-backup --diff
```

预览后追加 `--apply` 执行恢复。使用 TOML 时，把 `RESTORE_DIRECTORY` 放在
`[GRID_PARAMETERS]` 之前，先用 `RUN_MODE = "preview"`，再用 `RUN_MODE = "apply"`。
恢复模式忽略网格参数和备份目录。恢复完成后删除 `RESTORE_DIRECTORY`。

请用生成备份的 Suite 副本中的脚本执行恢复。脚本核对源码路径和校验和以保护后续改动。
完整修改计划在写入前生成，写入使用逐文件原子替换，遇到异常会回滚已完成的文件。
强制终止后可使用备份或保留的源码副本恢复。运行脚本时避免并发编辑同一源码树。

### 离线验证

```fish
../graspkit-tools/.venv/bin/python test/test_patch_grasp_grid.py
```

Fortran 编译器可用时，测试会编译最小用例，验证中央容量经转导后能用于共享数组。
