# F2018 标准兼容性首轮修复

2026-11-11：已否决的 strict/diag 独立求解器及其专用验收已[移除](rejected_rmcdhf_variants.md)。
下文涉及这些变体的修改、轮数和测试总数仅为历史记录，不属于当前支持范围。

2026-11-10 构建范围已收敛为普通 GRASP，CSFG 隔离；旧兼容 flags 已移除。
当前命令和 profile 以[构建策略](build_profiles.md)为准，以下数值记录保留原验证范围。

2026-11-09 更新：严格内部接口/NP、MPI Davidson 早收敛、普通 RANGULAR_MPI 空 rank
及四个 memory SETHAM 的 IOS 修复后，完成 strict 时序修复、32 配置矩阵及警告基线；
最新为[优先警告修复第一批](../reviews/f2018_warning_repairs_2026-11-09.md)。
默认 `f2018-audit` 保留静态存储；`f2018-check/release` 不保留，两者全量构建通过，
Linux 32/32 配置、累计 1,960 次 CTest 全通过；全开各 93/93。
3,046 文件基线及防新增检查通过；正确性复核 515→475，报告的 INTENT 不一致 9→0，完整 P0–P9 未完成。
strict 比较改为归一化后、阻尼前，原公式/阈值保持；1/2-rank 重启均 11 轮收敛。
同 O0/O2 的普通及 memory 一轮对照证据和此前 71/73 见
[Linux 历史报告](../reviews/f2018_linux_davidson_restart_2026-11-09.md)。此前 macOS 64/68 见
[上一轮报告](../reviews/f2018_strict_interfaces_npol_2026-11-09.md)。下面首轮结果属于 2026-11-08 的 Linux 基线。

日期：2026-11-08。修改前基线：`d1260ba404833d20d709b34d9f26b05c32a3b615`。
状态：当前 CMake 配置的全量目标已通过 GNU 15.3.0 的 F2018 编译和链接。
当前仅保留 `-fno-automatic`；后续接口批次已移除 `-fallow-argument-mismatch`，见
[接口验证记录](../reviews/argument_interfaces_2026-11-08.md)。本文其余批次表保留首轮实现记录。
数值与已知失败范围见[验证记录](../reviews/f2018_gfortran15_batches_2026-11-08.md)。

## 批次与保持的行为

| 批次 | 实施内容 | 保持范围 |
| --- | --- | --- |
| A，共享库与测试基础 | 普通/CSFG `nmtejj` 算术 IF 改为等值分支；CSFG `MIN(SIZE(...),n)` 统一 kind；重启计数及 QUAD 测试改为标准精度声明 | 表内容、匹配顺序、分配/复制顺序、64 位记录宽度与原测试容差 |
| B，平台服务 | 用标准 `iso_c_binding` 连接现有 POSIX 服务；保留 `sys_*` 入口，显式导入 SYSTEM/ACCESS/EXIT 的实现 | shell 返回状态、当前进程 cwd、成功/失败分支；新增对 mkdir 路径的正确引用 |
| C，计算程序标准语法 | 44 处 INT8 转为标准 INT；旧 kind 声明、维度声明顺序、已删除 DO、FORMAT、命令行 intrinsic 与内存统计转换 | 原计算次序、表达式、输入/记录格式和 solver 差异；不改变 SCF/RCI 算法 |
| D，验证与来源保护 | 新增五个单元/回归检查、固定一轮 SCF 对照脚本，以及可逆的来源适配清单 | 原始源码哈希和既有数值门槛不放宽 |

首轮语法批次未修改 `suite_parameters` 容量、网格初始化时序、EOL 权重、阻尼、停止规则、
NP 模型与物理常数。后续接口批次单独修正了 DDOT 返回类型、MPI handle 和数组实参；
本轮 strict SCNSTY 的归一化时序属于独立 F 类语义修复，保留原数值阈值。
大计数与其余内部接口问题仍需后续验收，不将这些问题混同于首轮语法等价变更。

## 平台服务接口

新模块 [suite_system](../../src/lib/lib9290/suite_system.f90) 归现有 `9290` 库所有，
现有 [sys_chdir](../../src/lib/mpi90/sys_chdir.f90)、[sys_getwd](../../src/lib/mpi90/sys_getwd.f90)、
[sys_mkdir](../../src/lib/mpi90/sys_mkdir.f90) 继续提供原调用接口。
这里使用的是标准 Fortran C 互操作加 POSIX 运行时，验证平台为 Linux；不宣称所有操作系统可直接使用。

- `run_shell` 保留 GNU SYSTEM 的原始 wait status，例如 shell `exit 7` 返回 1792。
  部分 CSFG 调用方通过 `USE suite_system, ONLY: SYSTEM => run_shell` 使用它；出现 SYSTEM 名称
  不再意味着隐式使用 GNU intrinsic。
- `change_directory` 改变当前进程目录；`working_directory` 获取并检查调用方缓冲区。
  现有调用契约仍为 0 成功、非零失败，不承诺不同运行时的具体 errno 数值相同。
- `path_access` 保留 POSIX 的存在性及 r/w/x 检查；调用方通过显式模块接口获得结果。
- `make_directory` 继续执行 `mkdir -p -m 775`，将路径作为单个经过引用的参数，支持空格、
  单引号和 shell 元字符。它不通过子进程 `cd` 模拟当前进程目录切换。
- `terminate_process` 保留原退出状态；省略参数为 0。测试覆盖非零退出与缓冲输出落盘。

## 标准 kind 与存储大小

计数转换采用现有 `vast_kind_param:LONG`，双精度采用 `DOUBLE`；明确四字节的索引声明
采用 `iso_fortran_env:INT32`。本机原有 32/64 位数据宽度保持，未使用全局整数扩宽选项。
将默认整数乘积最后转换为 LONG 并不能防止先前溢出；大计数边界仍需后续独立审计。

新模块 [suite_storage](../../src/lib/libmod/storage_bytes.f90) 使用 F2018 assumed-rank、
`STORAGE_SIZE` 和 `SIZE`，以 LONG 返回八位字节的 payload 大小。其调用对象须为有效的
已分配/已关联 intrinsic 标量或数组。CSFG 的 `report_mem` 继续只作内存诊断，原来注释掉的
调用点保持关闭；本次未启用该调试流程，也未修改其累计统计策略。

## 来源校验

原 [RMCDHF 来源清单](../reference/rmcdhf_upstream_snapshot.json) 保持不变。
[兼容性变更清单](../reference/rmcdhf_f2018_compatibility.json) 记录受保护文件中的
精确替换、次数和基线哈希；已扩展至 BLAS 声明、数组区段、MPI handle、内部接口和
SOLVE 参数意图、memory SETHAM 无效 IOS 判断移除，以及 GETALDWT 的输入长度修复，
逆向恢复后仍核对全部 492 个原始文件。
[check_upstream.py](../../test/npol/check_upstream.py) 先按行首精确匹配逆向恢复这些修改，
再剥离原有 NP hook，最后仍核对原始 hash。没有通过删除来源测试或覆盖历史 hash 放行。

负向检查已确认：额外修改阻尼表达式会改变原始 hash；撤回已登记的标准转换会使
适配清单匹配失败。后续变更必须另行记录其审核依据。

## 构建与验证入口

```sh
export FC=/usr/bin/gfortran-15
export OMPI_FC=/usr/bin/gfortran-15
cmake -S . -B build-gf15-f2018-batch \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
  -DMPI_Fortran_COMPILER=/usr/bin/mpifort \
  -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
  -DGRASP_BUILD_CSFG=ON -DGRASP_JJ2LSJ_OPENMP=ON \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DPython3_EXECUTABLE=../graspkit-tools/.venv/bin/python
cmake --build build-gf15-f2018-batch --parallel 1
ctest --test-dir build-gf15-f2018-batch --output-on-failure
```

2026-11-09 已隔离全部目标的模块输出目录，GNU 15.3 下全量 `--parallel 4` 构建通过；
上面的单任务命令保留为首轮复现入口。默认安装位置现为 `<build>/install`，支持显式前缀。
MPI 测试需允许本机进程通信。CMake 中独立编译的容量测试也启用 F2018，网格测试的
Fortran 编译器由配置中的 `CMAKE_Fortran_COMPILER` 传入。

固定一轮 SCF 的旧新对照入口见 [test/f2018/README.md](../../test/f2018/README.md)。
参考二进制必须来自独立保存的基线，同编译器和数学库配置；不要用当前源码重新构建后
覆盖旧基线。完整迁移阶段及后续验收仍以 [plan](../plans/f2018_migration_gfortran15_2026-11-08.md) 为准。

严格构建现在默认使用逐编译单元的 GNU JSON 诊断采集；已有无采集记录的构建需
新目录或 clean rebuild。`check.py` 核对源码/object/命令、编译器版本和配置独立的
固定诊断预算，新增/增量/缺记录均失败；不能以其他警告减少抵消。
见[逐文件索引](../reference/f2018/gfortran15-warning-baseline.md)。基线通过表示无超额，
`correctness_warnings_clear=false` 保留风险未清零状态，不将已登记警告视为安全豁免。

本批修复包括权重读取长度、预计算 Rk 积分统计、Davidson 接口与部分状态保留、
CSFG DACON 自有调试开关、INTERP 首点/短表/零值边界，以及被入口拒绝的 DGEMVMPI
非单位步长死代码。普通/memory 一轮对照和 INTERP 正常分支 71 点旧新输出保持一致。
边界缺陷用独立 oracle 验证；单点外推明确提示无法估计精度，不扩大 ACCY。
