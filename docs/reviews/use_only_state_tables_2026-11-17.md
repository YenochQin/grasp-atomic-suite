# USE ONLY 第六批：状态表与工作数组

日期：2026-11-17。基准为 HEAD `3e02e0f` 加此前未提交修改，接续
[状态模块与别名整理](use_only_state_aliases_2026-11-17.md)。本批保持接口、精度和计算行为。

## 范围

188 份源码补齐 280 条 ONLY，覆盖 35 个公开数据实体不超过六个、可在现有行宽内
完整列出的状态模块。以下数量较多的导入分别为：`tatb_c` 57、`pote_c` 32、
`stat_c` 26、`npot_c` 17、`mcp_c` 15、`ribojj_c` 12。

其余模块为 `facts_c`、`ribojj9_c`、`ribojj11_c`、`terms_c`、`ncc_c`、`l1_c`、`l2_c`、
`cnc_c`、`sbc_c`、`lic13_c`、`jqjc_c`、`sbdat1_c`、`corre_c`、`ovl_c`、`coun_c`、`left_c`、
`sacoef_c`、`cuto_c`、`titl_c`、`bcore_c`、`blim_c`、`keilst_c`、`stor_c`、`kkstart_c`、
`coeils_c`、`vinlst_c`、`vpilst_c`、`eigvec1_c`、`wchblk_c`。

ONLY 列出原有全部公开数据实体，保留 `double`、`byte`、`max_orbitals`、
`max_radial_points`、`radial_workspace_capacity` 等已有别名与容量常量。
不按局部使用次数删除导入，不改变量类型、维数、指针属性、初始化和数值表达式。
对编译模块导出表中模块自身、`iso_fortran_env`、`parameter_def` 的标识，确认实体种类为
MODULE 后才排除；它们重导出的参数仍保留。

范围仅限 `src/` 外部程序单元和单模块 CONTAINS 后的过程内 USE，模块头部重导出及
带重命名的导入保持不变。逐行检查只追加 ONLY，物理行数不变，没有新增超长行。

56 份导入 RMCDHF 文件追加
[state3 来源快照](../reference/rmcdhf_use_only_state3_snapshot.json)，绑定 state2 manifest；
不覆盖历史快照，来源检查继续验证前驱链与精确哈希。

## 验证

使用本机 Homebrew GNU Fortran 15.3.0，显式设置 FC/OMPI_FC；服务器工具链约定不变。
Python 使用工作区共享环境，采用独立 Release/Check 构建目录。

- 35 个模块的编译器导出表与完整 ONLY 列表一致；188 份文件的编辑范围与哈希核对通过。
- Fortitude C121 从 2,396 降至 2,116，默认诊断从 6,452 降至 6,172。
- 扩展正确性检查从 13,884 降至 13,604，差异恰为 280 条 C121；其余诊断按文件、规则、
  行号和消息比较完全一致。
- 来源链篡改负例通过，包含源码及历史 manifest 的改动。
- Release/Check 全量构建通过。208 个受影响 Release 编译对象逐字节一致，全部 1,204 份
  编译模块文件解压内容一致，接口未变。所有修改源码均列入当前构建目标，重复编译分别计数。
- 492 份导入 RMCDHF 文件的完整来源链及精确哈希检查通过。
- Release/Check 两套完整 CTest 均一次性 89/89 通过；`git diff --check` 通过。

本地忽略目录 `build-use-only-state3-20261117/` 保存修改前源码、导入/导出清单、哈希、
静态诊断、编译对象/模块接口对照、构建和回归日志。
没有新增物理计算；既有科学验证范围和 macOS memory 完整 SCF 限制仍见
[S3/S4 报告](static_numeric_interfaces_2026-11-17.md)。未提交、推送或安装覆盖生产二进制。

剩余 2,116 条 C121 继续按实际依赖、重命名、泛型和重导出契约逐批核定。
