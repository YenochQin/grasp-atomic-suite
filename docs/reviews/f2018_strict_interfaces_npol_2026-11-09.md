# F2018 严格接口、NP 路径修复与 GNU 15 验证

后续状态：本报告保留 macOS 上的历史结果；修复已提交为 `b8d6faf`。
Linux 后续修复与 71/73 验收、strict 未收敛基线见
[Linux 验证报告](f2018_linux_davidson_restart_2026-11-09.md)。

日期：2026-11-09。源码基线：`56769d91eae1e00889144e162458aceba0bdfe2a`，
结果对应其上的未提交工作树。本报告接续[验收设施报告](f2018_acceptance_facilities_2026-11-09.md)，
更新[迁移计划](../plans/f2018_migration_gfortran15_2026-11-08.md)的实施状态。

**严格 Check/Release 全量编译和链接通过；两套 CTest 各 64 通过、4 跳过，共 68 项。
强制验收仍返回失败，不能作为完整迁移或发布通过的结论。**
NP 在 audit/Check/Release 下各 15/15 通过；原版串行 RMCDHF 的同 O0、同 O2
一轮对照均保持七个二进制文件及全部科学文本/诊断一致。

## 工具链与范围

本机 macOS arm64，GNU Fortran 15.3.0、CMake 4.4.5、Open MPI 5.0.11、
OpenBLAS 0.3.34（LP64）。显式选择已安装的版本化 GNU 15，`FC/OMPI_FC` 指向同一编译器。
本机编译器路径与 Linux 参考路径不同；实际路径、编译器 SHA-256、MPI/数学库选择、
源码 hash、注册测试、CTest XML 和日志均保留在忽略提交的构建目录。

两个严格配置均打开 serial、MPI、CSFG、OpenMP；分别核对 **3,237 条编译命令**。
保留 `-std=f2018 -pedantic-errors -Wall -Wextra -Werror=implicit-interface
-Werror=implicit-procedure -fno-fast-math -ffp-contract=off`。
Check 为 O0、`-fcheck=all`、回溯及 invalid/zero/overflow 捕获；Release 为 O2。
没有使用参数放行、默认 kind 扩宽、全局 `-fno-automatic` 或删目标的方式通过严格构建。
audit 仍保留 `-fno-automatic`，用于回归 NP 原来的失败配置。

## 已实施修复

- **DINIT / MEM_FAIL 与内部接口**：调用方导入已有 `dinit_I`；`INIEST` 原先引用的
  `MEM_FAIL` 没有实现，现提供有显式签名的内部分配失败处理过程，输出对象、元素数和
  状态后 `ERROR STOP`。健康分配路径已有生产过程测试，未人为耗尽内存测试 OOM。
- **回调及公共库边界**：GDVD 使用显式 `dvd_operator` 回调契约；RKCO_GG、OUTSDA
  等回调使用匹配的过程接口。扩展公共 BLAS/LAPACK 接口的 SLAMCH、DLAMCH、LSAME、
  XERBLA；新增 `mpiu_I` 统一已有 MPI helpers 的声明，保留原通信实现。
- **应用和 CSFG 调用闭包**：复用已有 `_I`，缺失者按实际实现补齐并接入所属目标。
  修正 NCORE 的只读意图、MCPIN 的缺失参数、OVFLOW 的 SVEC 声明和 CSFG SETCSH
  名称；将非法的指针数组元素关联改为相应数组区段。程序内局部辅助例程放入 CONTAINS。
  修改范围较大是全量严格编译继续揭示后续调用闭包所致，未引入第二套科学库。
- **模块选择**：每个目标自己的模块目录始终置于依赖目录之前，包括已自行设置
  `Fortran_MODULE_DIRECTORY` 的目标，避免同名 `_I` 被其他应用的模块遮蔽。
- **Davidson 边界**：串行及 MPI 的 MGS_NRM 在最后一个新向量后不再为零列 DGEMV
  构造超出 B 上界的实参；零列运算本来没有数学效果。SCRA 的接口改为实际读写契约。
  本次没有更改 Davidson 收敛策略；额外的早收敛失败分支单列如下。
- **jj2lsj 初始化**：普通/CSFG 入口及系数测试在 SETISO 前调用 SETCON，先初始化
  AUMAMU 等单位常数，消除严格浮点检查发现的质量换算除零。
- **SOLVE 参数意图**：七个保留 solver 的 INV 改为 INOUT，JP/NNP 改为 OUT，
  实现与 `_I` 同步。INV 在调用者置零，某些分支保持原值；原 OUT 声明使该值在
  优化时失效。Release 对照直接捕获了 INV 输出从 0 变成 1，而七个二进制文件仍一致。
  修复后同 O2 的完整输出恢复一致；没有忽略该列或放宽文本比较。

容量 `215/590/600/127`、物理常数值、NP 模型、SCF 公式及停止条件保持原契约。
原版 RMCDHF 的接口修改以可逆差异登记到
[兼容变更清单](../reference/rmcdhf_f2018_compatibility.json)，原始来源 hash 未替换；
`check_upstream.py` **492/492 通过**。

## NP 摘要空路径

GNU 15.3 本机最小复现表明，`character(:), allocatable` 函数结果在
`-fno-automatic` 下传递的长度为零；先赋给 allocatable 临时变量也能复现。
同一小程序在 Check/O0、Release/O2 非静态配置下长度正常。这是本机实测，
不据此推断所有平台或 GCC 补丁版本的行为。

`npol_summary_path` 现使用依输入目录长度计算的显式长度字符结果，无固定路径上限，
保留原有斜杠/反斜杠及同目录输出语义。`npol_code_paths` 检查空输入、裸文件名、根目录、
空格、Windows 分隔符、300 字符目录、长短交替、直接实参和赋值后的精确长度/内容；
生产初始化另检查带空格子目录中的摘要落盘。

NP 入口检查更新为现行保留程序，并断言四个退役目录不存在；未删除四个原版入口的
NP hook 和顺序检查。原先因读取已退役目录失败的测试现已通过。

## 验证结果

| 检查 | 结果 | 范围 |
| --- | --- | --- |
| Check / Release 清单、编译命令、全量构建 | 两者通过 | serial/MPI/CSFG/OpenMP 全开，各 3,237 条编译命令 |
| Check 完整 CTest | 64 通过、0 失败、4 跳过 / 68 | 验收器退出 1，`strict_acceptance=false` |
| Release 完整 CTest | 64 通过、0 失败、4 跳过 / 68 | 验收器退出 1，`strict_acceptance=false` |
| NP audit / Check / Release | 各 15/15 通过 | 包含串行、真实 MPI、错误传播、路径及来源保护 |
| 生产 DDOT / CI 阻尼 | 每个严格配置 14/14 | 六个 solver；串行或 1/2/4 ranks，未冒充 CSFG DDOTMPI 覆盖 |
| Davidson 常规契约 | 每个严格配置 4/4 | 串行及 CSFG 库 MPI 1/2/4 ranks |
| jj2lsj | 每个严格配置 3/3 | 系数、排序和已有线程/回放回归 |
| RMCDHF 同 O0、同 O2 一轮对照 | 两者通过 | 每次七个二进制文件、全部科学文本/诊断精确一致 |
| 原版来源保护 | 492/492 | 逆向适配后仍匹配原始来源 |

新增 Davidson 契约调用真实 DINIT、INIEST、INIEST2、SPICMV2 和 GDVD。
检查步长与保护值、对角初猜、带非对角元的稀疏乘法，以及 6×6 三对角矩阵
（主对角 2、上下对角 -1）的最低本征值 `2-2*cos(pi/7)` 和残差，门槛 `1e-12`。
它不等于实际多根原子体系验证。SLAMCH/DLAMCH 的返回 ABI 另有机器常数检查。

四个跳过项为 `p0_trace_mpi_1/2`、`strict_mem_restart_mpi_1/2`，均为 Linux 专用路径。
未改 skip=77 或强制验收规则，未将平台不支持计作通过。

## 一轮旧新对照与保留的失败证据

目的为固定一轮软件等价性，预期同条件精确一致。输入为有限核 synthetic Z=10、
`1s^2`，`N=590, H=0.05, RNT=2e-7, HP=0, ACCY=1.5625e-8`，最大一轮，BLAS 单线程。
参考来自上述提交的独立 Git 归档，分别构建 O0 和 O2；O2 参考也显式使用
`-fno-fast-math -ffp-contract=off`。旧参考保留所需的静态存储，候选使用严格 profile。
本机 memory 版本依赖 Linux `/proc/meminfo`，因此显式选择 `--solvers rmcdhf`，
`full_serial_scope=false`。两种优化级别结果均达到预期，不表示 SCF 或物理模型已收敛。

首个 Release 对照误用 O0 参考，随后建立独立同 O2 参考；两次差异证据都保留。
同 O2 对照将差异定位到 SOLVE 的 INV 声明，修复后才得到通过结果。
二进制相等不能替代完整诊断比较。

另一个额外探针把 `[[2,1],[1,2]]` 嵌入 6×6 矩阵，其余对角为 4、5、6、7，
从第一坐标向量出发：串行通过，CSFG MPI 1/2/4 ranks 在严格配置下异常退出。
结合故障位置和源码，怀疑其 `NLOOPS > 2` 限制使极早的零残差继续进入 MGS 归一化。
这一原因尚需独立确认和修复，**本轮未解决该分支**；常规测试采用充分耦合的三对角矩阵，
明确不覆盖此失败分支。没有据常规 4/4 结果宣称整个 Davidson 求解器已验收。

## 复现与证据位置

从 suite 根目录执行；`GF15`、`MPI15_FC`、`MPI15_EXEC`、`OPENBLAS_LIB` 为本机已验证工具链，
共享 Python 使用 Tools 环境。Check/Release 分开构建，不能混用缓存：

```sh
export FC="$GF15" OMPI_FC="$GF15"
cmake -S . -B build-gf15-acceptance-check \
  -DCMAKE_Fortran_COMPILER="$GF15" \
  -DMPI_Fortran_COMPILER="$MPI15_FC" -DMPIEXEC_EXECUTABLE="$MPI15_EXEC" \
  -DBLAS_LIBRARIES="$OPENBLAS_LIB" -DLAPACK_LIBRARIES="$OPENBLAS_LIB" \
  -DPython3_EXECUTABLE="$(pwd)/../graspkit-tools/.venv/bin/python" \
  -DCMAKE_BUILD_TYPE=Debug -DGRASP_FORTRAN_PROFILE=f2018-check \
  -DGRASP_BUILD_SERIAL_APPS=ON -DGRASP_ENABLE_MPI=ON \
  -DGRASP_BUILD_CSFG=ON -DGRASP_JJ2LSJ_OPENMP=ON
../graspkit-tools/.venv/bin/python test/f2018/check.py build-gf15-acceptance-check --jobs 4
# 另一个目录使用 Release / f2018-release，同样运行 check.py。
../graspkit-tools/.venv/bin/python test/npol/check_upstream.py .
ctest --test-dir build-gf15-acceptance-audit -R '^npol_code_' --output-on-failure
../graspkit-tools/.venv/bin/python test/f2018/compare_first_sweep.py \
  --reference-bindir build-gf15-acceptance-baseline/build-release/bin \
  --candidate-bindir build-gf15-acceptance-release/bin \
  --output build-gf15-acceptance-release/audit/first-sweep-o2-inv-fixed \
  --solvers rmcdhf
```

对照输出目录须不存在；重复时另选新目录。正式强制验收在本机因四项 skip 退出 1。
MPI 回归需要本机进程通信权限。没有执行 install、提交或 SLURM 作业。

本地证据均在忽略提交的目录中：

- `build-gf15-acceptance-{check,release}/audit/acceptance/`：配置、构建、CTest XML、结果摘要和源码 hash。
- `build-gf15-acceptance-check/audit/first-sweep-inv-fixed/`：最终 O0 对照。
- `build-gf15-acceptance-release/audit/first-sweep-o2-inv-fixed/`：最终同 O2 对照。
- `build-gf15-acceptance-release/audit/first-sweep-current/`、`first-sweep-o2-matched/`：INV 修复前差异。
- `build-gf15-acceptance-audit/audit/npol-path-repro/`、`npol-current-test.log`：空路径最小复现和 audit NP 回归。
- `build-gf15-acceptance-check/audit/davidson-early-convergence-known-failure.log`：额外 MPI 探针失败。
- 两个严格构建的 `audit/before-entrypoint-jj2lsj-fix/`：NP 入口及 jj2lsj 修复前失败证据。

下一阶段仍需 Linux 上完成四项端到端测试和 memory 一轮对照、定位早收敛 MPI Davidson
失败、审计手写 INTENT 与跨调用状态、建立逐文件警告基线，以及补齐配置矩阵和 F00–F10。
当前构建仍有历史 INTENT 漂移、可能未初始化、静态数组边界等告警；编译成功不代表这些
风险已清零。`mpi_f08`、大计数、真实多根 SCF/RCI、跃迁、性能与物理验证均未完成。
