# RIS4 电子因子与 FICAL 两态同位素位移

2026-11-27，在 Suite `5a664b58` 加前批工作树修改的基础上，冻结了本批修改前的 `ris4`、`ris4_mpi`、`fical`。在 macOS arm64、GNU Fortran 16.2.0、Open MPI 5.0.11、OpenBLAS 0.3.34（LP64）上运行 Check 和 Release，各进程使用单线程 BLAS。目标是验证固定波函数的正常质量位移（NMS）、比质量位移（SMS）和场因子，以及这些量如何变成有符号的谱线位移。[生产入口回归](../../test/test_isotope_shifts.py)

## 单电子与多电子对照

氢样单电子夹具在 590 点网格上用动能平衡构造小分量。1s、2p₁/₂、2p₃/₂ 的领先 NMS 分别为 0.5000266253、0.1250005547、0.1250005547 a.u.，接近非相对论值 1/(2n²)。三个态的领先、相对论修正和总 SMS 均为零。接触密度依次为 0.3183095480、3.972741600e−7、1.212502936e−21 a.u.，分别对应 1/π、2p₁/₂ 小分量贡献和数值零。另用独立的 1s、2s 态文件验证交换输入文件时的上下态识别。[径向轨道夹具](../../test/test_radiative_transitions.py)

Ni I 固定 28 个电子、16 个轨道和 1,042 个 CSF，J=0、1、2、3、4 的偶宇称块分别保留 2、1、3、1、2 根。读取 MCDHF 混合流及同轨道 Coulomb+Breit CI 混合流，覆盖非零的二体反冲、接触密度、Seltzer 四个场因子和两个有效场因子。按报告的根、J 和宇称逐态读取。两种混合流各有 36 个不同态对。以下分项取自 MCDHF 中 J=0 第 2 根（上态）与 J=2 第 1 根（下态），能差 59,250.69 cm⁻¹。[Ni I 输入](../../data) · [RIS4 计算流程](../../src/apps/isotope/README.md)

| 谱线量 | 数值 | 单位 |
| --- | ---: | --- |
| ΔNMS | +204.480 | GHz u |
| ΔSMS | −137.744 | GHz u |
| ΔNMS + ΔSMS | +66.736 | GHz u |
| 质量位移 | −0.038354023 | GHz |
| 场位移 | −0.0058642614 | GHz |
| 总位移 | −0.044218284 | GHz |

表中的质量 60、58 u 和均方电荷半径差 +0.2 fm² 是指定的验证输入。电子因子按上态减下态，质量位移乘以 1/60 − 1/58。该线的两个有效场因子差为 −0.029323 GHz/fm² 和 +8.466e−6 GHz/fm⁴，分别乘半径差及其平方。球形 Fermi 模式另取两核 rms 电荷半径 3.80、3.77 fm，厚度均为 2.3 fm，以报告的二、四、六、八阶径向矩差重算场位移。测试同时核对二阶矩差与给定 rms 半径平方差。这里检验给定轨道、核模型和单位下的内部数值一致性，实验精度和活性空间收敛需要独立数据。

## 可复现的问题与修复

- **同尺寸错误基组。** 在仅有 1s、2s 的两 CSF 夹具中交换两行，让混合向量选中首行。旧程序复用旧顺序的角缓存，2s 接触密度得到 0.3183095480，冷算为 0.03978783520。现在一体和二体缓存同时核验完成记录、格式版本、基组摘要、算符种类、轨道标签基数及筛选阈值。交换基组后重新计算，再次复用时报告与重算逐字节相同。同一 Ni I 基组从接触密度模式切换到 Seltzer 模式可复用缓存，数值与冷算一致。[缓存身份](../../src/apps/isotope/open_isotope_shift_angular_cache.f90) · [缓存故障回归](../../test/ris4/check_cache_io.py)
- **输入顺序与场位移。** 将上表的同位素输入换为 58、60 u 及 −0.2 fm²，旧程序把该线的场位移算成 +0.0058642614 GHz。现在质量和场位移统一按重同位素减轻同位素报告，交换后该算例全部 36 行的质量、场及总位移与原顺序逐字节相同。接触密度模式仅含线性项，旧程序误用了相同系数作为二次项。整数形式的核电荷也已按实际值读入，原来的定宽格式将 `1` 读成了 `0.00001`。[FICAL 两态差](../../src/apps/isotope/isotope_shift_transition_main.f90)
- **报告不完整。** 旧程序接受缺少末行第二个有效场因子的报告。现在对能量、质量位移、密度和场因子各段检查记录、量纲和能级标识。回归包含截断、错误态号、错误场因子单位、混合 MHz u／GHz u 质量位移行，并检验失败时不生成新报告。[报告回归](../../test/test_isotope_shifts.py) · [质量位移单位回归](../../test/ris4/fical_units.py)

## 构建、MPI 与复验

每种构建各以串行及 MPI 1/2/4 ranks 运行四个固定输入（单电子、Ni I MCDHF 接触密度、Ni I MCDHF 拟合、Ni I CI 拟合）。32 份 RIS4 态报告均与各自基准逐字节相同。在相同构建内跨布局比较，全部 NMS、SMS 和接触密度在报告精度内一致。MCDHF 的 F₄ 最大绝对差为 1.0e−12 GHz/fm⁶，CI 的 F₆ 为 1.0e−14 GHz/fm⁸。对这 32 份输入的 FICAL 报告，按逐态因子独立重算并比较了 888 条谱线记录，即每种构建和布局各有 3 + 36 + 36 + 36 条。常规回归另覆盖交换同位素、两份态文件、球形 Fermi 核及 GHz、MHz、meV 三种输出单位。部署集群与形变核参数不在本次实测范围内。

Check 和 Release 各通过 14 项定向 CTest，另通过 5 项 F2018 数值、输入精度与验收驱动检查。从 Suite 根目录可复跑。工作区 `temp/tests/ris-p1/` 下的 `base/`、`candidate/` 保存输入、应答、日志和态及谱线报告，`faults/` 保存基准故障，`fixed/` 保存交换同位素后的复验，`comparison.json` 保存跨布局差值与表中算例。[P1 计划](../plans/plan-0004-f2018-follow-up-validation.md)

```fish
set -gx TMPDIR (realpath ../temp/tests)
for profile in check release
    ctest --test-dir "../temp/tests/rci-p1/$profile" --no-tests=error --output-on-failure \
        -R 'isotope_shifts|ris4_cache|fical_mixed|descriptor_dual_basis_pipeline|descriptor_v3_calculations|state_workflow|f2018_.*source_style'
end
```
