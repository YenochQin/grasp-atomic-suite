# 普通 GRASP 计算链导入与共享参数

日期：2026-11-08。代码接入基线：suite `179ca9ecfd96d54a5617d811b888299a611c13b3`。
状态：五类程序及共享容量接入已实施；这是 F2018 迁移前的准备，不是完整标准迁移或物理验证。

本页记录 d1260ba 导入快照；以下与上游逐字节相同的说明适用于导入时。
2026-11-11 普通程序目录已[去掉 90 命名](program_naming.md)；下表保留导入时的来源目录身份。
后续标准兼容性改动见 [F2018 修复说明](f2018_compatibility.md)；旧 orbopt 已[退役](orbopt_retirement.md)，
strict/diag 后续也已[移除](rejected_rmcdhf_variants.md)，当前 suite 为 19 个网格入口。

## 来源与范围

从工作区 `grasp/src/appl` 导入以下十个目录、815 个文件，保留原作者信息。
该来源副本没有可核验的独立 Git 元数据；[来源清单](../reference/standard_grasp_upstream_snapshot.json)
保存导入前各文件的 SHA-256 和路径。清单表示原始来源，不是适配后源码 hash。
没有复制 `.git`、计算输出、构建产物或第二套共享数值库。

| 源码目录 | CMake 目标 | 启用条件 |
| --- | --- | --- |
| `rangular90`、`rangular90_mpi` | `rangular`、`rangular_mpi` | 串行开关 / MPI |
| `rwfnestimate90` | `rwfnestimate` | 串行开关 |
| `rci90`、`rci90_mpi` | `rci`、`rci_mpi` | 串行开关 / MPI |
| `rbiotransform90`、`rbiotransform90_mpi` | `rbiotransform`、`rbiotransform_mpi` | 串行开关 / MPI |
| `rtransition90`、`rtransition90_phase`、`rtransition90_mpi` | `rtransition`、`rtransition_phase`、`rtransition_mpi` | 串行开关 / 串行开关 / MPI |

“串行开关”指 `GRASP_BUILD_SERIAL_APPS=ON`；MPI 目标要求 `GRASP_ENABLE_MPI=ON`
且找到 MPI Fortran。普通程序与已有 `*_csfg` 目标独立，不能互相当作别名。
`rhfszeeman95/rsms` 等不在本次五类程序的导入范围。

原应用 CMake 的源文件清单与算法保留，包括其中原先未启用的来源文件。
保留的旧 Makefile/BUILDCONF 用于来源追溯；suite 支持的构建入口仍是顶层 CMake。
所有目标复用现有 `mod/9290/dvd90/rang90/mcp90/mpiu90` 等所需库与 suite MPI helpers。

## 参数唯一来源

与 RMCDHF 相同，所有程序继续 `USE parameter_def`，由
[`parameter_def_M.f90`](../../src/lib/libmod/parameter_def_M.f90) 转导
[`suite_parameters_M.f90`](../../src/lib/libmod/suite_parameters_M.f90) 的六个常量：

| 参数 | suite 默认值 |
| --- | --- |
| KEYORB | 215 |
| NNNP | 590 |
| NNN1 | NNNP + 10 |
| NNNW | 127 |
| NNNWM1 / NNNWM2 | NNNW - 1 / NNNW - 2 |

没有增加第二个配置文件，也没有更改现有 suite 默认值。导入来源的 NNNW 为 100，
suite 统一使用已有的 127；文件内容和算法不因这个容量差异自动获得物理验证。

多数来源例程已经引用 `parameter_def`。本次仅在普通 `rwfnestimate90` 的六个
Fortran 文件进行容量适配：`frmtfp.f90`、`frmhyd.f90`、`solvh_I.f90`、`tail_I.f90`、
`frmrwf_I.f90`、`summry_I.f90` 移除局部 NNNP/NNN1 定义并改为模块导入。
其中接口的数组长度必须随中央配置变化，不能继续固定为 590/600。
其余导入文件与来源逐字节一致；没有改写数值求解、角系数、RCI 或跃迁算法。

运行时网格 `N/RNT/H/HP/ACCY` 保留各入口的上游初始化与交互顺序，不改为
统一执行的初始化例程。普通 RCI 保留原有公共网格初始化，不启用其注释中的点核分支。
编译期容量与运行时网格点数仍是不同参数；见[公共参数契约](common_parameters.md)。

## 网格补丁器接入

`atomic-suite` 布局新增八个网格入口：

- `rwfnestimate90/getinfo.f90`；
- `rci90/getcid.f90`、`rci90_mpi/getcid.f90`；
- `rbiotransform90/radpar.f90`、`rbiotransform90_mpi/radparmpi.f90`；
- `rtransition90/getosd.f90`、`rtransition90_phase/getosd.f90`、`rtransition90_mpi/getosdmpi.f90`。

另登记普通 RCI 两个 `lodres.f90` 续算 ACCY 位置。当前布局合计 25 个网格入口、
4 个 RCI 续算入口。`rangular` 没有径向网格初始化入口，但使用同一轨道容量/打包参数。
仅修改 `--nnnp` 仍只写中央配置文件；显式修改其他网格参数时，按原补丁器规则
更新对应入口，支持预览、应用、幂等检查、备份和恢复。

## 构建与测试接入

在已配置 gfortran 15 及匹配 MPI/BLAS 的环境中，从仓库根目录构建：

```sh
cmake -S . -B build-gf15-standard-import \
  -DCMAKE_Fortran_COMPILER=gfortran-15 \
  -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
  -DGRASP_BUILD_CSFG=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build build-gf15-standard-import --parallel 1
ctest --test-dir build-gf15-standard-import --output-on-failure
```

产物位于该 build 的 `bin/`。本次仍使用 suite 的 legacy 兼容选项，不是严格 F2018 构建。
MPI 包装器和库必须核实，不能只由 `gfortran-15` 的名称推断整个 MPI 栈的编译器版本。

现有 P0 trace/strict restart 小模型测试，在没有显式 `GRASP_PREPARATION_BINDIR`
且本次构建同时包含 `rwfnestimate` 和 `rangular_mpi` 时，自动使用本次 build 的准备工具。
首次配置即可注册，不需要先安装或等待可执行文件出现后再配置。
显式指定的外部准备目录仍按既有规则检查，不会自动换成别的版本。

新增 `standard_grasp_alternate_capacity`：在临时副本把中央 NNNP/NNNW 改为 790/149，
用 CMake 选中的 Fortran 编译器检查共享数组的 790/800 容量以及普通波函数接口的形状。
它使用过程指针关联验证真实接口；若 `tail_I` 仍硬编码 590，编译应失败。
生产参数文件不受该测试影响。

完整 RCI、双正交和跃迁数值回归仍按 [F2018 计划](../plans/f2018_migration_gfortran15_2026-11-08.md)
后续阶段实施。本次编译、容量和既有小模型回归不能证明实际原子体系的精度、活性空间收敛或性能。

## 2026-11-08 验证结论

GNU 15.3.0 完成全量构建和独立 serial-only 构建；35 项补丁器测试与非默认容量测试通过。
CTest 合并结果为 28/39 通过，11 项 NP/strict 重启失败已在未修改基线和原版准备工具上复现。
全量并行构建另触发既有 gj90 模块输出冲突，改用单任务构建完成。
详细范围、失败信息及对照方法见[验证记录](../reviews/standard_grasp_integration_validation_2026-11-08.md)。
