# F2018 验收设施、CI 阻尼分支与一轮对照

后续更新：本报告保留设施建立时的失败证据；DINIT/MEM_FAIL、NP 路径及后续回归问题
的修复结果见[严格接口与 NP 验证](f2018_strict_interfaces_npol_2026-11-09.md)。
当前 Check/Release 全量构建通过，各 64/68 测试通过、4 项 Linux 专用测试跳过；完整验收仍未通过。

日期：2026-11-09。来源基线：`56769d91eae1e00889144e162458aceba0bdfe2a`；
本报告对应其上的未提交工作树，不是发布版本。原 2026-11-08 Linux 报告保留为历史记录。

**结果：验收设施已实现，14 项生产 CI 阻尼测试通过；全量严格验收仍未通过。**
本轮未提交 SLURM 作业、未安装候选产物、未开展真实原子的收敛或性能研究。

## 构建与验收设施

- `cmake/FortranProfiles.cmake` 提供 `legacy`、默认 `f2018-audit`、`f2018-check`、
  `f2018-release`。均要求 GNU 15；Check 使用 O0/运行时检查/浮点异常捕获，Release 使用 O2，
  二者都启用标准及显式接口门槛，禁用全局隐式静态存储和快速浮点选项。
- BLAS/LAPACK 的 CMake 自带链接探针故意使用无参数隐式调用；只在包发现过程中暂不施加
  接口告警升级，随即恢复。suite 的实际编译命令仍逐条检查，生产目标没有该豁免。
- 每个目标独立输出 `.mod`；保留 jj2lsj 并发闭包及 `-frecursive`。
  `gj90/gj90_mpi` 不再在普通构建后复制到源码 `bin/`。显式安装前缀被保留，默认安装位置
  为 `<build>/install`，本轮没有运行 install。
- `configure.sh` 支持 profile、构建目录、安装前缀、串行及 MPI/CSFG/OpenMP 开关，
  使用 shell 数组保留参数中的空格，并拒绝未知参数。
- `test/f2018/required.json` 是**当前已实现支持范围**的目标/测试清单，不冒充 F00–F10。
  `check.py` 刷新配置后检查清单、实际编译器及逐条命令，再构建和运行完整 CTest。
  严格配置缺 MPI/Python/OpenMP 时失败；零测试、漏注册、未构建目标、失败和 skip 均不通过。
- 已注册 `f2018_unit/serial/mpi/properties/io` 标签。尚无真实体系夹具，故没有伪造
  `f2018_real_case` 测试。该标签及完整 F00–F10、警告基线和完整配置矩阵仍待后续验收。
- 本地 `audit/acceptance/` 保存配置、构建、CTest XML/日志、源文件 hash、编译器 hash、
  环境/库选择、测试与编译条目数量、失败状态。非严格 profile 不能获得严格验收结论。

## DDOT 消费路径及本轮发现的边界缺陷

`test/f2018/ci_damping.F90` 编译各应用的真实 `MATRIX/MATRIXmpi`，使用生产共享模块、
内存管理和实际 BLAS，仅控制矩阵装配、本征求解与内存角系数输入。没有复制阻尼算法。

夹具固定两个根、两个分量、非零块偏移及保护值，连续调用三次：

1. 首次调用得到单位基，确认 `FIRST` 绕过阻尼。
2. 下一次本征向量为 `(0.6,0.8)` 与 `(-0.8,0.6)`；第一根阻尼为 0、第二根为 0.5。
   第二根混合后与第一根有非零重叠，必经 DDOT/Schmidt 分支；解析预期经原符号约定后
   为 `(0.8,-0.6)`。检查分量、归一化、正交性、有限值及块边界，容差 `1e-12`。
3. 关闭阻尼，检查未混合的本征向量，作为控制路径。

两个串行版本及原版 MPI、mem MPI、strict MPI、diag MPI 的 1/2/4 ranks，合计 **14/14 通过**。
四 ranks 还覆盖无矩阵列的 rank。CSFG 使用独立 `DDOTMPI` 及展开布局，不在这 14 项中；
本夹具不是完整 SCF、多根物理模型或 Davidson 的端到端验证。

负向控制在隔离副本恢复旧 `INTEGER, EXTERNAL :: DDOT`，仅该错误副本允许隐式接口编译；
同一生产分支夹具在运行时按预期报告 `DDOT damping branch differs from analytic roots`。
正式目标没有兼容选项豁免。

开启运行时边界检查后，四个 MPI 实现先在 DDOT 之前暴露既有错误：`IENDC` 以 `1:NCF+1`
分配，却写入 `0:NCF`。本轮改为显式分配 `0:NCF`，不改变元素数量、矩阵公式或停止规则。
内存版本同时导入已有 `read03_mem_I`。原 RMCDHF 来源 hash 不改，新增可逆变更登记；
**492 个文件来源保护检查通过**。修复前失败日志和修复后结果分别保留。

## 一轮对照

对照目的为固定一轮计算的软件等价性。使用有限核 synthetic Z=10、`1s^2`、
`N=590, H=0.05, RNT=2e-7, HP=0, ACCY=1.5625e-8`，单线程 BLAS；
参考来自上述提交的独立 Git 归档，候选来自当前工作树，两侧 GNU 15.3/O0/OpenBLAS。
预期为二进制及保留的科学文本精确一致。

`compare_first_sweep.py` 现在要求恰好七个二进制输出、固定 MCP 记录形状、完整波函数/
混合系数记录、有限数值、程序完成标记及恰好一轮/迭代上限状态。分别保存并比较 stdout、
stderr、原生日志和 summary；只规范化明确列出的日期、时间及 Wall time 数值。
警告、能量、迭代与终止信息全部保留。失败也写 `comparison.json`，不依赖可被 `python -O`
移除的 assert。记录输入和可执行文件 hash，但不据此自动证明两套构建的全部环境一致。

实际 `rmcdhf` **7/7 二进制及全部科学文本/诊断一致**。macOS 上原版 `rmcdhf_mem` 在
读取 `/proc/meminfo` 时退出，故本轮实际对照显式选择 `--solvers rmcdhf`，报告为部分范围。
默认仍要求两个 solver；没有自动 skip 或用另一个 solver 替代。

脚本负向测试覆盖：两侧同时缺 MCP 文件/记录、空文件、截断记录、非有限诊断、未完成，
以及单侧警告、能量、迭代变化；均被拒绝。仅允许的时间元数据变化可以通过。

## 实测矩阵与未通过项

环境：macOS arm64，GNU Fortran 15.3.0，Open MPI 5.0.11，Homebrew OpenBLAS 0.3.34，
CMake 4.4.5，工作区共享 Python 3.15.2。按本机安装解析 GNU 15 的路径并显式设置
`FC/OMPI_FC`；没有复用旧 GNU 16 构建目录。MPI 初次受沙箱 socket 限制，随后在允许
本机进程通信的环境重跑；下表为重跑结果。

| 检查 | 结果 |
| --- | --- |
| audit，串行/MPI/CSFG/OpenMP 全开，`--parallel 4` | 全量编译链接成功；3147 个编译条目 |
| audit 完整 CTest | **50 通过、9 失败、4 跳过，共 63**；验收器返回失败 |
| legacy，串行开，MPI/CSFG/OpenMP 关 | 全量构建成功；CTest **21/26 通过、5 个 NP 失败** |
| 严格 Check 和 Release | 配置和命令检查通过；全量构建均被既有隐式接口挡住 |
| 严格 Check/O0 的串行阻尼两个目标及 BLAS 契约 | **3/3 通过**，不代表全套 strict 已通过 |
| 严格 Release/O2 的同三项 | **3/3 通过** |
| 配置负向检查 | 缺 MPI、缺 Python、参数放行、GNU 16 均拒绝 |
| Python 验收器/一轮对照脚本测试 | **7 项 unittest 通过**，每项含多个正反场景 |
| 来源保护 | **492/492 通过** |

九个 NP 失败仍为摘要空路径组；四个 skip 是 P0 trace 与 strict restart 的 1/2-rank
Linux-only 运行。严格构建首先在 `libdvd90/iniest.f90` 的 `mem_fail/DINIT`、
`lib9290/iniest2.f90` 的 `DINIT` 隐式接口处失败；本轮未扩大到全部内部接口修复。
这些与 2026-11-08 Linux 的 36/47 结果属于不同平台及测试集合，不直接比较通过率。

## 可追溯证据

以下均为忽略提交的本地目录；第三方需按文档重新生成，不把日志缺失视为测试通过。

- `build-gf15-acceptance-baseline/`：修改前提交的独立源码、四个串行二进制及初次输入/日志。
- `build-gf15-acceptance-audit/audit/acceptance/`：全开配置的当前构建与完整 CTest。
- 同目录下 `../before-iendc-fix/`：IENDC 越界反例；`../ddot-negative-control/`：旧 DDOT ABI 反例。
- `build-gf15-acceptance-audit/audit/first-sweep-rmcdhf-final/`：实际一轮对照的全部输入/输出。
- `build-gf15-acceptance-{check,release,serial}/`：严格/串行配置及其验收日志。
- `build-gf15-acceptance-config-negative/`：配置拒绝测试。

下一步优先修复并独立验收内部显式接口和 NP 摘要问题，在 Linux 恢复 memory/strict
完整运行，再补实际多根/RCI/跃迁与 F00–F10 基线。完整 P1/P3/P9 的退出条件保持不变。
