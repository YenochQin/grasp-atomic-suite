# USE ONLY 第五批：保留精度别名与容量常量的状态模块

日期：2026-11-17。基准为 HEAD `3e02e0f` 加此前未提交修改，接续
[小型共享状态模块整理](use_only_small_state_2026-11-17.md)。本批保持接口、数值精度和计算行为。

## 范围

216 份文件补齐 383 条 ONLY，涉及以下 15 个模块：

| 模块 | 数量 |
| --- | ---: |
| hblock_c | 72 |
| npar_c | 64 |
| default_c | 63 |
| fixd_c、mcpa_c | 各 24 |
| blkidx_c | 20 |
| damp_c、mcpb_c、peav_c、wfac_c | 各 16 |
| orba_c、lagr_c、invt_c、ncdist_c | 各 12 |
| node_c | 4 |

ONLY 显式列出各模块原有全部公开数据实体，包括状态变量、`double` 精度别名以及
`max_orbitals`、`max_radial_points` 等容量常量；不根据局部使用次数删减实体。
例如 `use npar_c` 改为：

```fortran
use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
```

通过修改前编译模块符号表核实 15 个模块的公开实体列表。对于表内模块自身及
`iso_fortran_env`、`parameter_def` 的元数据标识，确认其实体种类是 MODULE 后才从数据实体
列表中排除；不把 `double` 等来源于这些模块的参数误当作元数据排除。

范围仍是 `src/` 内外部程序单元及单模块 CONTAINS 后的过程内导入，保持模块头部重导出、
显式重命名导入和复杂接口不变。逐行核对只追加 ONLY，行数不变，没有新增超长行。
数值常量、变量声明、初始化值和调用顺序未改。

108 份导入 RMCDHF 文件追加
[state2 来源快照](../reference/rmcdhf_use_only_state2_snapshot.json)，绑定前一份状态模块
ONLY manifest；历史快照不覆盖，来源检查继续验证前驱链及精确哈希。

## 验证

使用本机 Homebrew GNU Fortran 15.3.0，显式设置 FC/OMPI_FC；服务器工具链约定不变。
Python 使用工作区共享环境，采用独立 Release/Check 构建目录。

- 15 份基准编译模块导出表与完整 ONLY 列表一致，公开变量和参数没有减少。
- Fortitude C121 从 2,779 降至 2,396，默认诊断从 6,835 降至 6,452。
- 扩展正确性检查从 14,267 降至 13,884，仅减少 383 条 C121；其余诊断按文件、规则、
  行号和消息比较完全一致。
- 216 份文件的编辑范围及修改后哈希核对通过；来源链篡改负例通过。
- Release/Check 全量构建通过。238 个受影响 Release 编译对象逐字节一致，全部 1,204 份
  编译模块文件解压内容一致。所有修改源码均列入当前构建目标，多目标重复编译分别计数。
- 492 份导入 RMCDHF 文件的来源链及精确哈希检查通过。
- Release/Check 两套完整 CTest 均一次性 89/89 通过；`git diff --check` 通过。

本地忽略目录 `build-use-only-state2-20261117/` 保存修改前源码、导入/导出清单、哈希、
静态诊断、编译对象和模块接口对照、构建/回归日志。
没有新增物理计算；既有科学验证范围和 macOS memory 完整 SCF 限制仍见
[S3/S4 报告](static_numeric_interfaces_2026-11-17.md)。未提交、推送或安装覆盖生产二进制。

剩余 2,396 条 C121 继续按实际符号来源、别名、泛型和重导出契约分批处理。
