# Landé 因子与混合态身份

2026-11-27，在 Suite `5a664b58` 及上一批跃迁修改的工作树上，分别保留修复前的 `gj`、`gj_mpi`，与修复后的同构建、同布局程序对照。目标是检验普通磁矩项、反常磁矩项及其和是否属于正确的原子态，并确认 MPI 分列累加后能恢复串行结果。环境为 macOS arm64、GNU Fortran 16.2.0、Open MPI 5.0.11、OpenBLAS 0.3.34（LP64）。Check 启用运行时和浮点异常检查，Release 使用 `-O2`，OpenMP 与 BLAS 均为单线程。[程序流程](../../src/apps/gj/README.md) · [生产入口回归](../../test/test_lande_factors.py)

## 单电子 LS 极限

取氢样 1s、2p₁/₂ 和 2p₃/₂ 大分量，用动能平衡构造小分量，在 590 点网格上联合归一化。纯态非相对论普通项分别趋于 2、2/3、4/3。对这组径向函数，联合归一化给出的普通项为 $g_{LS}/[1+(2nc)^{-2}]$，反常磁矩项为 $2\times0.001160\,(g_{LS}-1)$。两者独立计算期望值，因此可以检出角因子、修正符号或输出归一化的错误。[径向夹具](../../test/test_radiative_transitions.py) · [算符归一化](../../src/apps/gj/write_lande_g_factors.f90)

| 态 | 普通项 | 反常磁矩项 | 总 g 因子 |
| --- | ---: | ---: | ---: |
| 1s ²S₁/₂ | 1.9999733747 | +0.0023200000000 | 2.0022933747 |
| 2p ²P₁/₂ | 0.66666444787 | −0.00077333333333 | 0.66589111453 |
| 2p ²P₃/₂ | 1.3333288957 | +0.00077333333333 | 1.3341022291 |

## Ni I 九根的态匹配

固定仓库的 28 电子、16 轨道、1,042 个 CSF 和径向函数。J 从 0 到 4 的五个偶宇称块分别选择 2、1、3、1、2 根。描述符 SHA-256 为 `1a965eff0498dd3641f914f1fbdcb626fe8518045431c3ec6bb33b1b3b43acc5`。用这套轨道分别运行纯 Coulomb RCI 与 Coulomb+Breit RCI。后者开启横向光子项，频率缩放系数为 `1e-6`，真空极化、质量位移和自能均关闭。两次 CI 各自生成独立的 `.cimix`。[Ni I 基组](../../data/test_header.toml) · [MCDHF 参考报告](../../data/test.gj)

- 对每个 J、宇称块，枚举块内根的一一配对，取混合向量绝对重叠之和最大的一组。MCDHF 与纯 Coulomb CI 的最小绝对重叠为 0.9999999999999998。七个 J>0 态的总 g 因子在报告精度下一致，两个 J=0 态参与匹配并按程序约定略去输出。
- 加入 Breit 后，最小绝对重叠为 0.9996865652（J=2，第 2 根）。总 g 因子变化最大的是 J=2 第 3 根，从 1.4484219135 到 1.4555860710，普通项和修正项的变化分别为 +0.0071475763 和 +0.0000165812。这是同轨道下更换 CI 哈密顿量的结果。
- 反转块内根的存储顺序、略去 J=2 第 2 根并翻转各向量整体相位后，其余态仍按原根编号得到相同的三项。回归按报告精度采用绝对容差，普通项和总和各为 `5e-10`，较小的反常磁矩项为 `5e-12`。[逐态参考值与选根用例](../../test/test_lande_factors.py)

## 混合流身份错误

保持基组摘要、文件长度、偏移和系数正确，只把单电子 ²P₃/₂ 的块标签改为 J=5/2。修复前，串行 Check 正常结束，以错误的 J 归一化，给出总 g 因子 0.87337492136。错误宇称、错电子数、错块边界和重复根也能通过原读取器。现由 GJ 与 RHFS 共用的读取器对照已加载的基组，检查电子数、轨道数、总 CSF 数、块数、每块 CSF 数、J、宇称及根编号。串行和双 rank 的 GJ 回归逐一构造八类错误流，均非零退出并给出对应诊断。[读取器](../../src/apps/hyperfine/read_mixing_blocks.f90)

另从串行及双 rank RHFS 实际入口读入合法 Ni I 波函数，均生成对角及非对角报告。相同入口读入错标 J 的流时，均在矩阵元计算前退出。RHFS 会预先打开报告，这两次失败留下的 `.h` 和 `.hoffd` 均为零字节。超精细常数及带符号非对角元的数值对照归入下一批。

## MPI 归约与复验

对单电子、Ni I MCDHF、纯 Coulomb CI、Coulomb+Breit CI 四种输入，Check 和 Release 均运行串行及 MPI 1/2/4 ranks，32 份合法输入报告逐字节等于各自的修复前基准。三列数值在不同布局间的报告差均为零。

另在临时构建的串行及 MPI 2/4-rank 内核中记录普通、反常两个算符的对角约化矩阵元。24 次观测均复现未插桩报告。各 rank 的归约结果相同，用 `math.fsum` 独立重加局部矩阵元，与 MPI 归约的最大绝对差分别为 `6.94e-18`（MCDHF，J=4 第 1 根）和 `1.11e-16`（纯 Coulomb CI，J=3 第 1 根）。与串行矩阵元的最大绝对差分别为 `7.98e-17` 和 `4.44e-15`，均出现在纯 Coulomb CI 的 J=2 第 3 根。单电子模型只有三个 CSF，第四个 rank 的两项局部贡献逐态均为零。[MPI 分列累加与归约](../../src/apps/gj/compute_gj_only_mpi.f90)

本机记录位于工作区 `temp/tests/gj-p1/`。两组 RCI 应答、输入和输出在 `ni-input/`、`ni-breit/`，原始错标输入在 `baseline-wrong-J/`，32 组配对结果在 `comparison.json`，各 rank 的局部值及归约值在 `matrix/` 和 `matrix-summary.json`。以下命令从 Suite 根目录复跑已配置、已构建的两个 profile。各 12 项定向 CTest 全部通过。部署平台的完整构建矩阵及活性空间的物理收敛仍由后续批次验收。[实施计划](../plans/plan-0004-f2018-follow-up-validation.md)

```fish
set -gx TMPDIR (realpath ../temp/tests)
for profile in check release
    ctest --test-dir "../temp/tests/rci-p1/$profile" --no-tests=error --output-on-failure \
        -R 'lande_factors|descriptor_dual_basis_pipeline|descriptor_v3_calculations|state_workflow|f2018_.*source_style|f2018_mixing_output'
end
```
