# rmcdhf 局域核极化势：现行接口与代码验证

> 状态：代码已实施；2026-10-24 在 Linux/GNU Fortran 下完成下述代码测试。
> suite 基准：`7f9bf2d25ad1a2c6074c61c798cc686899677da4`；实现版本为首次引入本文的 Git 提交。
> RMCDHF 来源：2026-10-24 从工作区 `grasp/src/appl/rmcdhf90*` 新复制，逐文件 SHA-256 见文中清单。
> 验证范围：人工输入、数值内核、生产消费者、初始化与 MPI 协议。
> 未执行真实原子/离子 SCF、能级或同位素位移对照及物理收敛验证。

查询实现提交：`git log -1 --diff-filter=A --format=%H -- docs/implemented/rmcdhf_nuclear_polarization.md`。

## 支持范围

仅以下四个目标显式读取 NP 配置：

| 目标 | 源码入口 | 输入位置 |
| --- | --- | --- |
| `rmcdhf` | `rmcdhf/getscd.f90` | 当前输入目录的 `npol.inp` |
| `rmcdhf_mem` | `rmcdhf_mem/getscd.f90` | 当前输入目录的 `npol.inp` |
| `rmcdhf_mpi` | `rmcdhf_mpi/getscdmpi.f90` | `permdir/npol.inp`，仅 rank 0 读取 |
| `rmcdhf_mem_mpi` | `rmcdhf_mem_mpi/getscdmpi.f90` | `permdir/npol.inp`，仅 rank 0 读取 |

已有 `rmcdhf_orbopt*`、独立 strict/diagnostic 求解器、性质程序和外部 RCI 不读取该配置。
公共 `NUCPOT` 不读取文件；库被链接不等于应用启用了 NP。
`isodata`、交互 stdin 顺序、`rwfn.out` 和 `rmix.out` 格式保持原状。

## 原版来源与保留边界

原版四个应用复制到 `src/appl/rmcdhf/`、`rmcdhf_mem/`、`rmcdhf_mpi/`、
`rmcdhf_mem_mpi/`；原有 `src/appl/rmcdhf90*` 继续构建 `rmcdhf_orbopt*`。
此前误加到 orbopt 的 16 个文件改动已撤回，并逐字节核对为 suite 基准版本。
NP 的生产消费者测试现在编译新复制的原版应用例程。

[来源清单](../reference/rmcdhf_upstream_snapshot.json) 记录 492 个复制文件的原路径与
SHA-256，其中 480 个为 Fortran 源码。源 `grasp/` 当前没有独立 Git 元数据，不能把
父工作区的 Git SHA 当成其上游版本，因此以本次实际读取的文件快照作为可核查基线。
`npol_code_upstream_copy` 剥离明确列举的 NP 接入修改后，核对全部 492 个文件。
修改仅在每个变体的 `GETSCD*`、接口、主调用点及 `DACON` 的显式接口导入中。
`SCF/SCFMPI`、`IMPROV*`、`SOLVE`、正交化、节点判据和矩阵例程保留复制源内容，
不引入 orbopt 的额外控制与轨道更新流程。

公共库继续使用 suite 的单份实现。与本地 grasp 的比较显示，既有差别是集中参数
转导、`COUNT` 的可选 `GRASP_TRACE_COUNT` 日志及其上下文变量、MPI 帮助例程扩展
和注释；节点计数数值部分、`NUCPOT`、`RINTI` 与 `libdvd90` 相同。
这些已有共享扩展没有在此次迁移中修改；没有复制第二套公共库。
这一区分意味着应用求解源码来源于原版，不宣称整个 suite 可执行文件与外部 GRASP
二进制逐字节相同，也没有进行两者的真实原子结果或耗时比较。

## `npol.inp` v1

文件缺失时关闭 NP，不提问、不生成 NP 输出文件。显式关闭可写：

```fortran
&npol
  schema_version = 1
  enabled = .false.
/
```

文件存在时必须是一个以 `/` 结束的 `&npol` namelist，允许空白、`!` 注释及
Fortran 引号转义；拒绝额外块/内容、未知字段、未知版本、错误类型与无法读取的文件。
文件大小上限为 65536 字节。字符串先读入能容纳完整文件的缓冲区，再检查长度，
避免 namelist 的静默截断。

| 字段 | 要求 |
| --- | --- |
| `schema_version` | 文件存在时显式为 `1` |
| `enabled` | 默认 `.false.`；显式 `.true.` 才启用 |
| `use_e1`, `use_e2` | 默认 `.false.`；启用 NP 时至少选一项 |
| `units` | 启用时为 `atomic`；最多 16 字符 |
| `model_id` | 启用时为 `flambaum2021_scalar_e1_e2`；最多 64 字符 |
| `parameter_source` | 启用时非空，最多 1024 个 Fortran 字符（UTF-8 按字节计），不含控制字符；记录来源和 E2 约定 |
| `isotope_z`, `isotope_a` | 启用时为显式整数；`Z>0` 且与运行核电荷相同，`A>=Z` |
| `k1`, `b1` | E1 开启时必填：有限的 `k1>=0` 和 `b1>0` |
| `k2`, `b2` | E2 开启时必填：有限的 `k2>=0` 和 `b2>0` |

负哨兵区分未提供参数与显式零强度。未选择的通道贡献严格为零，不使用它的
数值参数；整体关闭时只检查语法、版本和字符串表示，不使用通道参数。
`A>=Z` 是输入一致性检查，不是核素存在性验证；程序不从平均原子量推断质量数。
内部不联网、不拟合参数、不换算 fm/MeV/B(EL)。输入生产者负责物理参数与约定。

以下**仅为软件夹具**，不表示任何核素的物理参数：

```fortran
&npol
  schema_version = 1
  enabled = .true.
  use_e1 = .true.
  use_e2 = .true.
  units = 'atomic'
  model_id = 'flambaum2021_scalar_e1_e2'
  parameter_source = 'artificial software fixture; E2 convention: test only'
  isotope_z = 2
  isotope_a = 4
  k1 = 16
  b1 = 2
  k2 = 128
  b2 = 2
/
```

## 势表示和生命周期

实现的模型为

\[
V_1=-\frac{K_1}{r^4+b_1^4},\qquad
V_2=-\frac{K_2}{r^6+b_2^6},\qquad
ZZ=ZZ_{\mathrm{Coul}}-r(V_1+V_2).
\]

`r,b1,b2` 使用 `a0`，`k1` 使用 `Eh a0^4`，`k2` 使用 `Eh a0^6`。
输入强度已包含模型所需系数，内核不再乘 `e²α/2`。
E2 极化率与截断参数必须来自同一约定。
参数确定、跃迁方向与单位换算见[参数准备参考](../reference/nuclear_polarization_parameters.md)，
其中区分论文固定 50 keV 的 E2 修正极化率、实验跃迁能及截断拟合的适用条件。

调用顺序是 `RADGRD → NUCPOT → NPOL_INIT[_MPI] → SETRWF*`。
核电荷、Fermi 参数、`NNUC`、网格、`ACCY` 时序不改变。

- [npol_C](../../src/lib/libmod/npol_C.f90) 保存受保护的配置、Coulomb 快照、E1/E2/总势和状态。
- [npol_input](../../src/lib/lib9290/npol_input.f90) 提供 `NPOL_READ_CONFIG` 和 `NPOL_VALIDATE_CONFIG`，返回状态，不终止进程。
- [npol_kernel](../../src/lib/lib9290/npol_kernel.f90) 提供 `NPOL_BUILD_GRID`、`NPOL_COMPOSE`；只写活动区 `1:N`，失败不提交输出数组。
- [npol_initialization](../../src/lib/lib9290/npol_init.f90) 编排串行初始化和摘要输出。
- [npol_mpi](../../src/lib/mpi90/npol_init_mpi.f90) 提供 MPI 初始化、显式字段广播和非零退出适配器。

初始化先 reset，再捕获调用者刚经 `NUCPOT` 生成的纯 Coulomb 势；重复合成始终
从保存的快照开始。重新初始化前，调用者必须重新生成 Coulomb 势，不能把旧总势当基线。
失败保留调用者的本轮 Coulomb 数组，`npol_ready=.false.`，NP 配置关闭；成功后配置
及数组在 SCF 内不再更新。关闭/全零强度走直接复制路径。

内核通过二进制尾数、指数和缩放半径求值，避免先计算 `r**4`/`r**6` 所导致的
中间溢出；检测不可表示的势与合成溢出，保留可表示的次正规数。
原点不用 `ZZ/r`，负半径、非有限数、非法 `N` 或短数组均返回错误。

`YPOT`、`DACON` 和 `RINTI(...,0)` 继续读取同一份完整 `ZZ`；`RINTI(...,1)`
仍是动能路径。没有第二份 NP 矩阵项；冻结轨道、EOL 权重、角向系数、二电子项和
SCF 停止条件的实现未修改。四份 `DACON` 仅补充了已有 `DPBDT` 的显式接口导入。

## MPI 与输出

所有 rank 先集体检查活动长度，再由 root 解析、校验和构建势；广播配置中的每个字段
及完整 E1/E2/总势数组。各 rank 使用本地纯 Coulomb 快照合成，所有合成成功后才由
root 输出摘要，再发布各 rank 的状态。错误消息来自首个失败 rank 并广播给所有 rank；
适配层返回错误，应用调用 `MPI_Abort(...,1)`。不使用 root 单独调用 `stopmpi` 的路径。

每个 rank 存储完整 `ZZ`。原有 MPI `YPOT` 的 `ZZ/NPROCS` 与后续求和归约保持原样。
新路径使用可分配字符串，不把拼接后的文件名塞进固定长缓冲区；当既有 `permdir`
达到其 128 字符容量时明确失败，避免沿用可能已截断的路径。既有目录帮助例程未改写。

启用时在配置文件同目录写 `npol.summary`，含版本、开关、核标识、来源、单位、
参数、`N` 和 `composed = T`。未选择通道的未提供值以内部负哨兵记录，应结合开关解读。
打开、写入、截断、关闭及写后内容校验均检查状态；MPI 同步输出失败。
关闭时不创建或改写摘要；原有旧摘要不被自动清理，不能用其存在判断本次是否启用 NP。

## 构建与验证记录

2026-10-24 使用 GNU Fortran 16.2.1、CMake 4.4.4、Open MPI 5.0.11，
保持默认容量 `NNNP=590`。两套 Debug 构建均完成，MPI 构建确认四个新支持目标及四个旧 orbopt 目标均存在。
NP 生产例程与测试驱动启用了 GNU 数组边界检查。没有执行安装或修改共享 Python 环境。

从仓库根目录运行：

```bash
cmake -S . -B build-npol-upstream-serial -DCMAKE_BUILD_TYPE=Debug \
  -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-npol-upstream-serial -j4
ctest --test-dir build-npol-upstream-serial -L '^npol_code$' --output-on-failure

cmake -S . -B build-npol-upstream-mpi -DCMAKE_BUILD_TYPE=Debug \
  -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-npol-upstream-mpi -j4
ctest --test-dir build-npol-upstream-mpi -L '^npol_code$' --output-on-failure
```

| 验证集 | 结果 |
| --- | --- |
| 串行 `npol_code` | 9/9 通过 |
| MPI 配置 `npol_code` | 14/14 通过，含真实 1/2-rank 启动 |
| 串行既有纯代码回归 | 7/7 通过 |
| MPI 配置既有纯代码回归 | 10/10 通过 |

既有回归选择了 `lib9290_quad`、`shared_grid_capacity`、`strict_scf_logic`、
`rmcdhf_pair_transaction_logic`、`rdensity_natural_orbitals`、
`rmcdhf_orbopt_round_state_logic`、`grasp_grid_patcher`；MPI 配置另含
`mpi90_local_sparse`、`start_trace_parity`、`p0_trace_projection`。
未运行全量 CTest，未运行依赖原子数据的回归。

本机沙箱禁止 MPI 所需的本地 socket；最终 MPI 测试在允许本地通信的执行环境中运行。
MPI 启动失败和超时不会被当作预期错误测试通过。
构建保留原有代码的 legacy FORMAT/旧 MPI 参数类型警告；没有将这些警告视为物理验证。
本次全量 MPI `-j4` 构建曾触发既有 RHFS 串行/MPI 目标对同名 `.mod` 文件的并发写入冲突；
移除这次失败留下的两个 `getmixblock_I.f90.o` 构建对象后，使用
`cmake --build build-npol-upstream-mpi -j1` 完成构建。该性质程序的构建竞态未在本次修改中修复。

本轮筛选测试日志保存在各新构建目录的 `npol-code-ctest.log` 和 `related-code-ctest.log`。
此前 `build-npol-serial` / `build-npol-mpi` 中的旧产物不作为本轮验证或运行入口。

详细测试结构见 [test/npol/README.md](../../test/npol/README.md)，设计依据见
[原实施计划](../plans/rmcdhf_nuclear_polarization_implementation.md)。
这些结果支持软件契约与接入正确性，不支持物理精度结论。下一阶段需要另行设计并授权
目标核素参数、SCF/轨道响应、目标可观测量与收敛验证。
