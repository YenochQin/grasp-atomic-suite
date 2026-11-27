# 超精细常数与非对角矩阵元的相位

2026-11-27，`rhfs` 和 `rhfs_mpi` 在本机通过 A、B 及带符号非对角矩阵元验收。基准是 Suite `5a664b58` 加上前批已验证的工作树修改，冻结了改动核矩换算前的两个入口。Check 和 Release 均使用 macOS arm64、GNU Fortran 16.2.0、Open MPI 5.0.11 及 OpenBLAS 0.3.34（LP64）。Check 启用运行时与浮点异常检查，Release 使用 `-O2`，每个进程的 BLAS 与 OpenMP 各用一个线程。[生产入口回归](../../test/test_hyperfine_constants.py)

## 单电子解析对照

取氢样 1s、2p₁/₂ 和 2p₃/₂，沿用 590 点网格及动能平衡小分量。用单电子径向积分和角因子独立计算 MHz 常数，用总角动量 F 的耦合系数计算 Hartree 矩阵元。在 I=1、磁矩 μ=1 核磁子、电四极矩 Q=1 barn 时，得到以下结果。[径向夹具](../../test/test_radiative_transitions.py) · [单位常数](../../src/lib/core/physical_constants.f90)

| 电子态或有序态对 | A（MHz） | B（MHz） |
| --- | ---: | ---: |
| 1s ²S₁/₂ | 254.42504235 | 0 |
| 2p ²P₁/₂ | 10.601149278 | 0 |
| 2p ²P₃/₂ | 2.1202298557 | 3.9160796395 |
| ⟨2p ²P₃/₂，2p ²P₁/₂⟩ | 0.66257182990 | 0.84783918175 |

在夹具的正系数相位下，最后一行的两个 F 分量分别为 −9.0970325247e−10 Hartree（F=1/2）和 +5.5787892701e−10 Hartree（F=3/2）。将核矩分别改为 (μ,Q)=(−2,0)、(0,3)，单独检验磁偶极和电四极对这两个分量的贡献。另以 I=1/2、μ=1、Q=0 核对 A 的反比缩放。所有对角 F 分量均从报告中的 A、B 经 Casimir 公式独立重算，换算系数为 6.5796839207e9 MHz/Hartree。反宇称态对的相互作用为零。

## Ni I 九根、根顺序和相位

固定库内 Ni I 的 28 个电子、16 个轨道、1,042 个 CSF，偶宇称的 J=0、1、2、3、4 块各取 2、1、3、1、2 根。使用夹具指定的 I=1、μ=1 核磁子、Q=1 barn。所测量的是这组给定核矩和电子波函数的性质。分别读取 MCDHF `.mix`、同轨道纯 Coulomb 和 Coulomb+Breit RCI 的 `.cimix`。Breit 算例将横向光子频率缩放为 `1e-6`，关闭真空极化、质量位移和自能。[核参数](../../data/isodata) · [CI 选项与描述符摘要](evidence-0029-lande-g-factors.md)

- **MCDHF 与纯 Coulomb CI**。每个 J、宇称块内以绝对向量重叠匹配根，再将每个矩阵元乘以两端向量重叠的符号。最小绝对重叠为 0.9999999999999998。七个 J>0 态的对角 A、B 均复现仓库参考值。计入非对角常数，A、B 最大绝对差分别为 2.0e−9、4.0e−9 MHz。报告中共有 76 个 F 分辨矩阵元，其中 55 个为非对角元，最大绝对差为 1.0e−18 Hartree。27 组非对角态对覆盖 ΔJ=0、1、2，同 J 态对另核对了两个方向的厄米对称性。[对角参考](../../data/test.h) · [态对与 F 分量参考](../../data/test.hoffd)
- **改变存储顺序和相位**。把真实 Coulomb CI 的各块根行逆序，只对 J=2 第 3 根的整个向量变号。它与 J=1 第 1 根之间 F=1 的矩阵元由 −3.3078567821e−8 变为 +3.3078567821e−8 Hartree。串行和双 rank 的报告均在相位对齐后恢复原值。回归比较所有态与态对，使用绝对容差 2e−7 MHz、2e−17 Hartree 及相对容差 2e−10，包含接近零的矩阵元。[匹配及符号回归](../../test/test_hyperfine_constants.py)
- **加入 Breit**。九根仍可一一匹配，最小绝对重叠为 0.9996865652。最大的对角 A、B 变化都在 J=2 第 3 根，A 从 467.74594505 降至 460.26281355 MHz，B 从 631.87644198 升至 634.30455538 MHz。F 分辨矩阵元的最大变化为 5.3810e−9 Hartree，来自 J=2 第 2 根与 J=1 第 1 根的 F=2 分量。这些差值对应哈密顿量的改变。

## 零自旋与工作区

对合法的 I=0、μ=Q=0 输入，原程序在 A 的 MHz 换算及磁偶极约化核矩中均除以 I。Check 的串行及 MPI 1/2/4 ranks 异常退出，Release 返回成功但写出 NaN。共用的换算过程现将零自旋核的超精细作用置零。两种构建的四种布局均给出有限的零 A、B 和空的 F 矩阵元表，电子 Landé 因子仍可读。负自旋输入单独给出诊断。串行例程现在也在返回前释放两块按能级对分配的 Landé 工作数组。[核矩换算](../../src/apps/rhfs/hyperfine_nuclear_factors.f90) · [串行工作区](../../src/apps/rhfs/calculate_hyperfine_matrix_elements.f90)

## MPI 归约与复验

对单电子、Ni I MCDHF 和两种固定轨道 CI，在 Check、Release 上各运行串行及 MPI 1/2/4 ranks，共 32 组非零自旋计算。每组的对角和态对报告均与修复前同构建、同布局的报告逐字节相同。跨布局的 A、B 和 F 分量在报告精度内一致。

在另外 24 次观测中，串行及 MPI 2/4 ranks 的临时二进制额外记录了归约前后的五路电子约化矩阵元。用 `math.fsum` 重加局部贡献，与 rank 0 归约结果的最大绝对差如下。与串行最远的一项为 Ni I MCDHF 的 J=4 第 2 根对角电四极元，绝对差 1.24e−14。单电子算例分配到四个 rank 时，最后一个 rank 的所有局部贡献为零。插桩前后的全部报告逐字节相同。[MPI 累加与归约](../../src/apps/rhfs/calculate_hyperfine_matrix_elements_mpi.f90)

| 归约分量 | 最大绝对差（归一化前） |
| --- | ---: |
| 磁偶极 ΔJ=0 | 2.78e−17 |
| 磁偶极 ΔJ=1 | 6.94e−18 |
| 电四极 ΔJ=0、1、2 | 各 4.44e−16 |

工作区 `temp/tests/rhfs-p1/` 下的 `base/` 和 `candidate/` 保存两版的输入、日志和报告，`phase/` 保存真实 RCI 及改换根顺序和相位后的混合流，`matrices/` 保存各 rank 的原始矩阵元。`comparison.json` 和 `matrix-summary.json` 分别汇总报告与归约差异。以下命令从 Suite 根目录复跑已构建配置。Check、Release 的 14 项定向 CTest 各自全部通过。部署平台与其他活性空间仍需独立验证。[P1 实施进度](../plans/plan-0004-f2018-follow-up-validation.md)

```fish
set -gx TMPDIR (realpath ../temp/tests)
for profile in check release
    ctest --test-dir "../temp/tests/rci-p1/$profile" --no-tests=error --output-on-failure \
        -R 'hyperfine_constants|lande_factors|descriptor_dual_basis_pipeline|descriptor_v3_calculations|state_workflow|f2018_.*source_style|f2018_mixing_output'
end
```
