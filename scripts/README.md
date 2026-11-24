# GRASP 辅助脚本

2026-11-16：网格工具同时识别原版及语义变量名，命令行/TOML 参数保持不变。
`atomic-suite` 现覆盖十五个普通网格入口和两个 RCI 续算入口，CSFG 专用入口已移除。
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

## 保存 RMCDHF 计算结果

使用固定文件名的 RMCDHF 计算结束后，可在启动计算的目录运行 [save_rmcdhf](save_rmcdhf)。
按[状态名计算](../README_ZH.md#输入与输出)已直接生成后续程序的输入。以下归档示例假设 suite 的 `bin/` 在 `PATH` 中。也可直接运行源码目录中的脚本。脚本使用系统 `sh` 和 GNU coreutils。

```fish
save_rmcdhf results/Zn_odd_AS1
```

目录的最后一段是状态名，以下用 `name` 表示。脚本创建独立副本，保留工作文件。若目标已存在或复制失败，命令以非零状态退出，已有归档保留。

| 计算目录中的文件 | 归档目录中的文件 | 内容 |
| --- | --- | --- |
| `rcsf_desc.parquet`、`rcsf_header.toml` | `name_desc.parquet`、`name_header.toml` | 原生 CSF 描述符与头文件按原字节成对复制，保留 header 哈希绑定。 |
| `rwfn.out` | `name.wfn` | 径向轨道，保留 Fortran 非格式化记录 |
| `rmix.out` | `name.mix` | OL/EOL 的 CSF 混合系数，保留 Fortran 非格式化记录 |
| `rmcdhf.sum` | `name.rmcdhf.summary` | 计算摘要 |
| `rmcdhf.log` | `name.rmcdhf.log` | RMCDHF 交互记录 |
| `rangular.log` | `name.rangular.log` | 角系数程序的交互记录，存在时复制 |
| `isodata`、`npol.inp` | `isodata`、`npol.inp` | 存在时按原名复制的核参数和核极化输入 |

脚本依据摘要中的计算模式处理混合系数。OL/EOL 缺少 `rmix.out` 会报错，(E)AL 会跳过目录中遗留的 `rmix.out`。归档轨道文件为 `name.wfn`，OL/EOL 混合系数文件为 `name.mix`。有 `isodata` 的 OL/EOL 归档可整体移到另一目录，再从归档中运行 `gj`，输入同一状态名并选择 MCDHF 混合系数。供 `gj` 等程序使用时，状态名宜用 1 至 23 个无空格的 ASCII 字符。

MPI 通过 `disks` 单独配置永久 I/O 目录时，从含有 `rmcdhf.log` 的启动目录执行以下命令。脚本从 `--io-dir` 读取描述符、轨道、混合系数、摘要和核数据。

```fish
save_rmcdhf --io-dir /path/to/serial-io results/Zn_odd_AS1
```

## 导入原版 GRASP 的 rsave 结果

在原版 `rsave` 结果所在目录，运行 [import_rsave](import_rsave) 并传入当时使用的状态名。
命令安装到 Suite 的 `bin/`，使用 Python 3.10 或以上版本和支持 V3 的 `rcsfs`。

```fish
cd /path/to/grasp-results
import_rsave Zn_odd_AS1
```

`rcsfs` 默认从 `PATH` 查找，也可用 `--rcsfs /path/to/rcsfs` 指定。

| 原版结果 | 同目录生成的 Suite 文件 |
| --- | --- |
| `name.c` | `name_desc.parquet`、`name_header.toml`，保留 CSF 行顺序和对称性块 |
| `name.w` | `name.wfn` |
| `name.m`、`name.cm` | `name.mix`、`name.cimix`，分别在源文件存在时生成 |
| `name.sum`、`name.log`、`name.alog` | `name.rmcdhf.summary`、`name.rmcdhf.log`、`name.rangular.log`，分别在源文件存在时生成 |

脚本要求同一状态的 `.c` 和 `.w`。摘要标明 `(E)AL` 时跳过原目录中遗留的 `.m`。
轨道和混合系数原样复制，要求原版输出与 Suite 的 Fortran 二进制记录格式兼容。已有 `isodata`、LSJ 标签和能级表继续在原目录使用。
原始文件与同名 Suite 文件受保护，目标文件已存在时命令报错退出。

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

脚本检查十五个网格入口、两个 RCI 续算入口和六项共享容量的转导，并扫描整个 `src/`，
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
