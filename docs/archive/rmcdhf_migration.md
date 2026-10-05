# rmcdhf_test 合入 grasp-atomic-suite

> 文档状态：已实施迁移的历史记录、旧默认值已被替代。
> 核查日期：2026-10-05；源码基准：`8b298a6`。
> 四个 orbopt 目标和共享库已落地；正文中性质程序 590 点默认值、多入口补丁及旧测试数量属于迁移阶段。当前默认 N=NNNP 与集中脚本布局见[公共参数](../implemented/common_parameters.md)和[脚本说明](../../scripts/README.md)。
> 全部文档状态见[分类索引](../README.md)。

本文记录初次迁移的取舍。后续已将公共参数集中管理，所有有限核程序默认
点数统一为 `NNNP`；当前设计与验证以 [common_parameters.md](../implemented/common_parameters.md)
为准。下文的 590 点默认值和脚本布局描述属于迁移当时的状态。

轨道优化开发现在由 `grasp-atomic-suite` 维护，与性质计算程序共用一套数值库。
原 `rmcdhf_test` 工作目录和 Git 历史保留作为参考。

## 来源与迁移范围

- suite 基线：`e6ca561ded6b0c1c98e228194eeb24e532a28a27`。
- 轨道优化来源：`rmcdhf_test` 的
  `38f6bba8dcaf735e08bc64f9dcdd1362caf24747`，加上迁移时工作区中的
  `orbopt` 重命名、网格修改工具及其测试。
- 迁入全部四个 `src/appl/rmcdhf90*` 目录（498 个源代码/构建文件）、
  已跟踪的 `test/` 文件、网格修改脚本和测试；未复制 `.git`、编译产物或大型计算数据。
- 合并采用源码迁移，两个远程仓库的提交历史和远程设置保持独立。
  新代码应在 suite 提交；没有自动提交、推送或删除原仓库。

## 公共库与数值约定

`libmod`、`lib9290`、`libdvd90`、`mpi90` 只有一份，轨道优化与性质程序共用。
保留 suite 的 `mpiu.f90`，包括 `mpix_startup`、`mpix_shutdown` 和 root reduction
辅助函数；轨道优化继续使用其中原有的 `startmpi2/stopmpi2`。
`iniestmpi`、`spicmvmpi` 的执行逻辑一致，保留 suite 版本。
加入 `COUNT` 的可选诊断和 `coun_C` 的轨道/进程标识；节点计算公式未改变。

全局容量统一为 `NNNP=2990`、`NNN1=3000`。保留原来的默认实际网格：

| 程序 | 有限核默认实际点数 |
| --- | --- |
| 四个 `rmcdhf_orbopt*` | `N=NNNP`，当前为 2990 |
| RHFS、RIS、g_J | `N=MIN(590,NNNP)`，当前为 590 |

g_J 使用 RHFS 的 `gethfd.f90`；轨道文件加载及交互输入仍可以改变实际网格。
这次迁移不将已有性质计算自动改成 2990 点计算。
外部 GRASP 的准备/RCI 工具仍需另行安装，并支持计算所用网格容量。

四个轨道优化命令为 `rmcdhf_orbopt`、`rmcdhf_orbopt_mem`、
`rmcdhf_orbopt_mpi`、`rmcdhf_orbopt_mem_mpi`；输入输出约定沿用 RMCDHF。
串行目标由 `GRASP_BUILD_SERIAL_APPS` 控制，MPI 目标由
`GRASP_ENABLE_MPI` 和可用的 MPI Fortran 工具链控制。
统一构建使用 CMake；迁入的 legacy Make/BUILDCONF 文件保留，不新增 suite 的 legacy 总构建入口。

## 网格修改工具

`scripts/patch_grasp_grid.py --layout atomic-suite --grasp .` 覆盖四个轨道入口、
RHFS/RIS 入口以及全局容量；g_J 自动使用同一 RHFS 入口。
默认 `grasp2018` 布局继续支持完整原版 GRASP。
两种布局都要求完整的已知入口，并扫描所有 Fortran 源文件；缺少文件或未知初始化
会使计划失败。默认预览，应用时备份，恢复时检查校验和及后续修改。
详见 [脚本说明](../../scripts/README.md)。

## 本次验证

- 开启 MPI 和全部串行程序：10 个可执行文件均成功构建并安装。
- CTest：公共库积分、共享网格容量、MPI 稀疏缓冲区、轨道事务、轮次状态、
  网格修改器共 6 项通过；网格修改器包含 12 个离线用例。
- 独立关闭 MPI 的构建成功，5 项 CTest 全部通过。
- 同一 `data/test` 输入在 suite 基线与合并版本的串行和双进程 MPI g_J 计算中，
  `test.gj` 输出逐字节一致；合并版本串行与 MPI 输出也逐字节一致。
- 真实合并 Fortran 源码副本上，网格脚本的应用、检查、重复应用和恢复通过；
  恢复后所有源码字节与原始副本一致。完整原版布局预览也通过。
- 迁入的 shell 脚本通过语法检查；新增 Python 修改通过 Ruff 检查。

没有运行大型 RMCDHF 原子计算或集群作业；这些仍需要仓库外的输入、外部 GRASP
准备程序及目标集群环境。历史结果中的具体编译目录/程序名保留作为记录，
机器专用绝对路径在迁入副本中用变量替代。

验证构建命令（先配置本机 Fortran/MPI/BLAS 工具链）：

```sh
cmake -S . -B build-merge-check -DCMAKE_BUILD_TYPE=Debug \
  -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-merge-check --parallel 8
ctest --test-dir build-merge-check --output-on-failure
cmake --install build-merge-check

cmake -S . -B build-merge-serial-check -DCMAKE_BUILD_TYPE=Debug \
  -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-merge-serial-check --parallel 8
ctest --test-dir build-merge-serial-check --output-on-failure
```

本次发现 149 个旧 `.mod` 文件残留在 suite 的 `src/` 下。它们会被 GNU Fortran
优先读取，造成新参数配合旧数组容量。已将这些生成文件移动到本地忽略目录
`build-merge-check/legacy-module-backup/` 后全量重新编译；CMake 增加了拒绝源码树
残留 `.mod` 的检查，CTest 也验证不同模块实际使用的容量一致。
