# 密度矩阵、自然轨道与性质输入身份

2026-11-27，以 Suite `bd4dcaa2cae8578ff315c3f8c76cdd098e1396f4` 为基准，在 Linux x86-64、GNU Fortran 16.2.1、LP64 OpenBLAS 0.3.34 和 Open MPI 5.0.11 上验证本批工作树。源码及测试补丁的 SHA-256 为 `e7a36089c326a906c38c70e058bf528bd8e1ca47652157639c41ecf10a943feb`。基准和候选均启用串行与 MPI、关闭 OpenMP，使用 `f2018-check` 和单线程 BLAS。候选另做相同布局的 Release 定向回归。[实施计划](../plans/plan-0004-f2018-follow-up-validation.md)

## 密度与自然轨道

三组冻结输入分别是带 2:4 统计权重的单电子模型、具有两重占据简并的单电子模型，以及仓库中 28 电子、16 轨道、1,042 个 CSF、五个偶宇称块的 Ni I 九根混合态。完整检查覆盖矩阵对称性、迹与占据总和、κ 壳层占据上限、本征残差和旋转矩阵正交性。用独立读取器核对 `RAD--WFN` 中的轨道身份、网格、P/Q 和原点系数。单电子还用自然轨道及占据重建原始径向密度。[生产回归](../../test/rdensity/regression.py) · [Ni I 基组](../../data/test_header.toml)

- 简并夹具用两个正交根张成与 `(1,1,1)` 垂直的平面。解析占据为 `0, 1/2, 1/2`，用该平面的投影矩阵比较简并子空间。候选在串行及 MPI 1/2/4 ranks 的最大投影元素误差为 `2.22e-16`，读回轨道的径向 Gram 矩阵相对单位阵最大误差为 `9.02e-16`。四 ranks 时第四个 rank 无 CSF 列，仍参加矩阵归约。
- 候选的单机串行和 MPI 1/2/4 ranks 共计算 12 组，每组的 `.d` 和 `.nw` 共 24 份文件均与各自的基准同布局结果逐字节相同。Ni I 不同布局间的矩阵和占据数最大绝对差各为 `1.78e-15`，读回后的自然轨道 P/Q 最大差均为零。这里使用固定波函数，跨节点和其他活性空间另行验收。
- MPI 根进程用 `MPI_IN_PLACE` 直接归约二维矩阵。各 rank 仍以宽整数计算元素数再检查窄化。自然轨道求解前检查对称性，并在分配之前验证 `DSYEV` 返回的工作区大小。原有排序、位相和返回状态继续由单元回归检查。[MPI 归约](../../src/apps/rdensity/density_parallel_mpi.f90) · [自然轨道求解](../../src/apps/rdensity/natural_orbitals.f90)

## RIS4 和 FICAL

只把单电子混合流目录中的 `3/2−` 改为 `5/2−`，保持描述符摘要、能量和混合系数。基准 RIS4 和 FICAL 生成了标成 `5/2−` 的结果，GJ 拒绝同一输入。本批由 RIS4 和 RDENSITY 共用基组检查，读取系数前逐项核对头记录的维度、块边界、J 和宇称。串行及 MPI 负例在生成计算报告前退出。[混合目录校验](../../src/lib/core/mixing_file_format.f90) · [RIS4 负例](../../test/test_isotope_shifts.py)

FICAL 现在按 `isodata` 的质量数选择参考核。保留正的 Fermi 占位参数、将质量数置为零时，旧报告写出 Fermi 核及 `2.085824 fm` 的参考均方根半径。新回归要求报告点核并把参考半高半径、均方根半径和皮厚写为零。[FICAL 读取](../../src/apps/isotope/isotope_shift_transition_main.f90) · [点核回归](../../test/ris4/fical_units.py)

首轮比较器仅规范化两种主机空闲内存行中的实时数值。把原 Linux 日志中的 `6.53` 和 `6.51 Gb` 两份输出重新送入规范化函数后，比较一致，原始日志和七份二进制摘要可继续追溯。[首轮比较器](../../test/f2018/compare_first_sweep.py)

## 复验

工作区 `temp/tests/den.lCOuzuM7/` 存有基准与候选的构建日志、`candidate.patch`、`frozen/` 输入、`density-runs-final/` 输出及 `density-summary.json`。`base-wrong-j/`、`check-wrong-j/`、`base-point-grid/` 和 `check-point-grid/` 保存两个性质输入问题的新旧实测。

Check 的定向清单共 19 项，包含密度、同位素位移、Landé、RHFS、RCI 续算、首轮比较驱动和源码检查。首跑通过 17 项，串行与单 rank 续算因缺少夹具生成程序 `rwfnestimate` 跳过，补建后两项均通过，详见 `check-tests.log` 和 `check-rci-retry.log`。Release 的 11 项定向测试全部通过，见 `release-tests.log`。在 Suite 根目录对已构建的配置运行以下命令，可重跑两个 profile 共用的密度和同位素位移检查。

```fish
set density_run (realpath ../temp/tests/den.lCOuzuM7)
set -gx TMPDIR "$density_run"
for profile in check release
    ctest --test-dir "$density_run/$profile" --no-tests=error --output-on-failure \
        -R '^(rdensity_|isotope_shifts_|fical_mixed|f2018_first_sweep_driver)'
end
```
