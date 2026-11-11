# F2018 兼容性分批修复验证

> 2026-11-11 状态更新：文中 strict/diag 独立求解器已按用户决定[移除](../implemented/rejected_rmcdhf_variants.md)；相关目标、命令和测试数量仅作对应版本的历史记录。

日期：2026-11-08。修复前基线：`d1260ba404833d20d709b34d9f26b05c32a3b615`。
本轮在其上保留先前的 `-std=f2018` 开关，实施[标准兼容性修复](../implemented/f2018_compatibility.md)。
本报告记录各批次提交前的工作树验证，不代表已经发布的版本。

**结论：当前 CMake 配置的全量 F2018 编译/链接通过；CTest 33/44 通过，余下 11 项与
已确认的旧基线失败集合完全相同。五项新增回归全部通过；四个串行 solver 的一轮 SCF
共 28 个文件与旧版本逐字节一致。**

## 1. 构建配置和范围

- GNU Fortran 15.3.0；`FC=/usr/bin/gfortran-15`，`OMPI_FC=/usr/bin/gfortran-15`。
- CMake 4.4.5，系统 OpenBLAS，Open MPI 5.0.11。
- Debug、普通串行程序、MPI、CSFG、已有 OpenMP 均开启；全量构建使用单任务。
- 实际选项含 `-std=f2018 -fno-automatic -fallow-argument-mismatch -g`；jj2lsj 保留原独立存储/线程配置。
- 所有 3588 个 CMake 编译条目均指向 gfortran 15 并带 F2018；候选构建没有使用 legacy 构建的模块搜索路径。
- 未安装新工具链、未调用 install、未提交 SLURM。既有 gj90 的 post-build 复制行为保留；
  单独保存的 `build-gf15-standard-import` 基线没有被重编译或覆盖。

全量候选在新目录 `build-gf15-f2018-batch` 从 F2018 失败状态分批重建，最终
`cmake --build build-gf15-f2018-batch --parallel 1` 返回 0。
另在全新 `build-gf15-f2018-serial` 关闭 MPI/CSFG/OpenMP，六个普通串行目标和两个新增
单元目标编译成功，两个相应单元测试通过。

这一结果不等于移除参数兼容选项后的严格标准认证；MPI/BLAS 旧接口仍会产生已知警告。
也没有把未导入的 `rhfszeeman95/rsms` 等计划目标计入“当前全量目标”。

## 2. 分批结果

| 批次 | 变更和复现 | 验证结果 |
| --- | --- | --- |
| A | 原共享库/QUAD 在 F2018 下先复现失败；改算术 IF、CSFG kind 和测试精度声明 | 普通 NMTEJJ、CSFG reallocation、QUAD 三项先通过；CSFG NMTEJJ 在相关库恢复后也通过 |
| B | 先前真实对象的 CHDIR/GETCWD 链接失败及 SYSTEM 编译失败 | 平台模块和 sys_* 真实接口通过；shell wait status、路径/错误、进程退出/文件缓冲测试通过 |
| C | 解开依赖后继续全量编译，逐项清理 INT8、旧 kind、维度顺序、DO、FORMAT、SIZEOF、命令行和 EXIT 扩展 | 当前配置全量编译和链接成功；保留静态状态和参数放行策略 |
| D | 原源码校验先拒绝合法标准转换；增加可逆适配清单并重新核对旧 hash | 492 个受保护导入文件通过；额外阻尼修改和撤回兼容转换的负向检查被拒绝 |

没有改变哈密顿矩阵公式、SCF/Davidson 更新流程、EOL 权重、轨道冻结集合、默认容量、
网格/ACCY 时序或原数值验收阈值。文件 I/O 类型宽度仍为本机原来的 32/64 位。

## 3. 新增与保留的回归

### NMTEJJ 精确对照

参考结果由 d1260ba 的普通/CSFG legacy 库、GNU 15.3.0 生成；两份结果一致。
固定 [reference](../../test/f2018/nmtejj_reference.json) 保存源 hash、输出 hash 与 22,701 行数。
输入覆盖 J=1..37 的奇数分支、J=9 两种 NK/ND 情况、表项命中/未命中及 J>37 上界返回。
候选普通和 CSFG 实现均与该结果逐字节一致；不是仅比较最终能量。

### 分配、尺寸与路径

- LONG 尺寸的整数/实数组扩容、缩容保留原前缀值；重启计数仍为 64 位。
- 存储字节检查涵盖 packed、默认整数、LONG、双精度、标量、空数组和非连续切片。
  `storage_bytes` 只对有效的已分配/已关联数据定义接口；未启用原来关闭的 REPORT_MEM 调用。
- 路径检查使用真实 `sys_*` 和新平台服务，覆盖空格、引号、shell 元字符、已存在目录、
  不存在目录、穿过普通文件的 mkdir、短 cwd 缓冲区及失败后 cwd 保持。
- GNU SYSTEM 的 `exit 7` 返回 1792 的原始 wait status 保留；进程 EXIT 的状态 0/7 和缓冲写出均通过。
- 既有 QUAD 解析积分与自然轨道、密度、jj2lsj、CSFG 文件/容量等检查继续使用原门槛。
  网格补丁器仍通过其 35 个用例；编译器由 CMake 配置传入，避免误用系统的另一个 gfortran 版本。

### 固定一轮 SCF 的旧新对照

目的：验证这批语法/接口适配保持固定迭代的计算行为。输入为软件测试使用的有限核
Z=10、`1s^2` 单 CSF 模型，固定 `N=590, H=0.05, RNT=2e-7, HP=0, ACCY=1.5625e-8`。
参考和候选使用相同 gfortran 15、数学库和单线程设置，分别运行 TF 初猜、串行角系数和一轮 SCF。
预期为二进制输出一致；不以“跑完一轮”判断收敛。

| solver | 比较产物 | 实际结果 |
| --- | --- | --- |
| `rmcdhf` | `rwfn.inp/rwfn.out/rmix.out` 与四个 MCP 文件 | 7/7 逐字节一致 |
| `rmcdhf_mem` | 同上 | 7/7 逐字节一致 |
| `rmcdhf_orbopt` | 同上 | 7/7 逐字节一致 |
| `rmcdhf_orbopt_mem` | 同上 | 7/7 逐字节一致 |

入口为 [compare_first_sweep.py](../../test/f2018/compare_first_sweep.py)，完整命令见
[测试说明](../../test/f2018/README.md)。结果满足本轮等价改写的预期，但不覆盖 CI 阻尼、
多根/近简并 SCF、完整 RCI、跃迁参数或实际原子的活性空间收敛。

## 4. 全量 CTest 的 11 项既有失败

在允许本机 MPI socket 的环境运行完整 CTest，44 项中 33 项通过、11 项失败，
没有跳过项。失败名集合与[导入阶段独立旧基线复现](standard_grasp_integration_validation_2026-11-08.md)完全一致：

- `strict_mem_restart_mpi_1/2`：既有小模型在要求的重启流程中未产生 `STRICT_CONVERGED`。
- NP 共九项：compose、output、entrypoints、四个 consumers、1/2-rank broadcast；
  同样因摘要文件路径为空而失败。

本批没有改变 NP 初始化或 strict 停止算法来使测试变绿；这些问题仍需独立诊断。
新平台路径测试不等于修复 NP 摘要路径为空的问题，也不把 strict 的未收敛状态解释为物理结论。
本轮运行结果必须表述为“没有新增回归失败”，不能表述为“全套数值验收通过”。

## 5. 后续工作和证据位置

下一批优先建立明确的 BLAS/LAPACK 接口、核对 DDOT 返回类型和 MPI datatype/count，
再逐步撤除 `-fallow-argument-mismatch`。每个可能改变结果的真缺陷使用独立参考验证。
SAVE/内存所有权、完整编译配置矩阵及真实多根模型按 plan 后续阶段推进。

本地证据位于忽略提交的 `build-gf15-f2018-batch/audit/`：分批 build/test 日志、
`ctest-full.log`、`ctest-summary.json`、`first-sweep-parity-with-logs/comparison.json` 及全部对照输入/输出。
最终一轮对照另核对 stdout 只有一次外层 SCF iteration，独立保存 stdin/stdout，保留程序原生日志。
serial-only 配置证据在 `build-gf15-f2018-serial/audit/`。
原失败诊断 `build-gf15-f2018-probe/audit/` 和 legacy 基线均保留，不以新结果覆盖旧记录。
