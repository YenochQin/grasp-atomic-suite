# 局域核极化势参数：文献定义与原子单位换算

> 核查日期：2026-10-24；代码基准：`5aa6ad8b1b6c88c42edf52ebbdc8be4cc05b4e20`。
> 本文是参数准备参考；未运行原子 SCF，未建立核素参数库或物理精度结论。
> 现行输入契约见 [RMCDHF NP 接口](../implemented/rmcdhf_nuclear_polarization.md)。

## 模型与原始来源

采用 Flambaum、Samsonov、Tran Tan、Viatkina，*Nuclear polarization effects in
atoms and ions*, Phys. Rev. A **103**, 032811 (2021)：
[期刊记录](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.103.032811)、
[作者预印本 PDF](https://arxiv.org/pdf/2102.00355)、
[预印本 HTML](https://arxiv.org/html/2102.00355v1)。
公式号以下按该论文；预印本和正式版公式号相同，表号不同：

| 内容 | arXiv v1 | 正式版 |
| --- | --- | --- |
| E1 数值截断半径 | IX | VII |
| E1 截断拟合系数 | V | VIII |
| E2 数值截断半径 | X | IX |
| E2 截断拟合系数 | VI | X |

正式版表号也与[作者上传的正式版全文](https://www.researchgate.net/publication/348958162_Nuclear_polarization_effects_in_atoms_and_ions)核对。
本轮 APS PDF 端点不可访问，公式同时以 arXiv PDF/HTML 核查。

代码存放的是完整系数：

\[
V_1=-\frac{K_1}{r^4+b_1^4},\qquad
V_2=-\frac{K_2}{r^6+b_2^6},\qquad
K_1=\frac{e^2\alpha_0^{E1}}2,\quad
K_2=\frac{e^2\bar\alpha_0^{E2}}2.
\]

这里 `b1` 对应论文 \(b\)，`b2` 对应 \(\tilde b\)。
它们通过 Eq. (42) 的能移匹配确定；不是核半径。
现有代码提供一套所有轨道共用的局域势，不自动选择电子参考态；若输入 \(ns\) 列参数，
就是采用论文讨论的 s 势近似，应在 `parameter_source` 中记录。
目前没有各 \(\kappa\) 通道专用的截断参数。
[依据：Sec. IV，Eqs. (40), (42), (48), (52)](https://arxiv.org/html/2102.00355v1#S4)。

## 强度：先明确跃迁方向

论文统一用退激发强度 \(B_\downarrow(EL)=B(EL;L\to0)\)。
NNDC 页面通常给 \(B_\uparrow(E2)=B(E2;0^+\to2^+)\)，且单位为
\(e^2\,\mathrm{barn}^2\)。由详细平衡：

\[
B_\downarrow(E1)=B_\uparrow(E1)/3,\qquad
B_\downarrow(E2)=B_\uparrow(E2)/5.
\]

这是角动量简并度转换，不是拟合修正。
[NNDC 定义及单位](https://www.nndc.bnl.gov/be2/help.html)。

### E1：巨偶极共振近似

论文 Eqs. (24), (25), (48) 给出

\[
E_{\rm GR}=95(1-A^{-1/3})A^{-1/3}\ \mathrm{MeV},
\]

\[
B_\downarrow(E1)=\frac{3}{8\pi}
\frac{Z(A-Z)e^2}{AE_{\rm GR}m_p},\qquad
\alpha_0^{E1}=\frac{8\pi}{3}\frac{B_\downarrow(E1)}{E_{\rm GR}}.
\]

上式的 \(\hbar=c=1\) 不应直接与 fm 混用。恢复单位后，定义
\(D_1=B_\downarrow(E1)/(e^2\mathrm{fm}^2)\)，可直接算

\[
D_1=\frac{3}{8\pi}\frac{Z(A-Z)}A
\frac{(\hbar c)^2}{(m_pc^2)E_{\rm GR}},
\]

其中两种能量用 MeV、\(\hbar c\) 用 MeV fm；右端数值即 \(D_1\)。
这是 TRK 集中在一个共振上的模型估计。若有实验强度，可记录其来源并替代强度；
若同时更改激发能或强度谱，现有截断拟合的匹配条件也需重新检查。
[依据：Sec. III B](https://arxiv.org/html/2102.00355v1#S3.SS2)。

### E2：与截断拟合配套的修正极化率

论文 Eq. (53) 使用固定参考能量：

\[
\bar\alpha_0^{E2}=\frac{8\pi}{5}
\frac{B_\downarrow(E2)}{\bar E_{\rm rot}},\qquad
\bar E_{\rm rot}=50\ \mathrm{keV}=0.050\ \mathrm{MeV}.
\]

它与 Eq. (56) 的 \(\tilde b\) 配套；不能保持这个 \(b_2\) 却把分母任意换成
实验 \(E(2_1^+)\)。若改用真实静态极化率，需要按 Eq. (42)/(55) 重新匹配截断半径。
[依据：Sec. IV D，Eqs. (52)–(56)](https://arxiv.org/html/2102.00355v1#S4.SS4)。

优先采用评价的实验 \(B(E2)\)。若只有集体形变参数，论文 Eq. (33) 与 NNDC
定义给出同一个换算：

\[
D_2\equiv\frac{B_\downarrow(E2)}{e^2\mathrm{fm}^4}
=\frac15\left(\frac{3Z}{4\pi}\right)^2
R_0[\mathrm{fm}]^4\beta_2^2,
\quad R_0=1.2A^{1/3}\ \mathrm{fm}.
\]

若采用 NNDC 的 \(S=B_\uparrow(E2)/(e^2\mathrm{barn}^2)\)，则
\(D_2=10^4S/5\)。因为 \(1\ \mathrm{barn}=100\ \mathrm{fm}^2\)，平方单位转换为
\(10^4\)，不能只乘 100。
[NNDC 形变定义](https://www.nndc.bnl.gov/be2/help.html)。

## 换成 `npol.inp` 的原子单位

为避免在 Gaussian/SI/自然单位之间重复乘电荷，最直接的流程是先把跃迁强度与
能量转换成原子单位，再计算 \(K\)。设

\[
a=52917.7210544,\quad h=2.7211386245981\times10^{-5},
\]

分别为 \(a_0\) 的 fm 数值与 \(E_h\) 的 MeV 数值；
\(m_p/m_e=1836.152673426\)，\(m_pc^2=938.27208943\) MeV，
\(\hbar c=197.3269804\ldots\) MeV fm。
[NIST 2022 CODATA](https://physics.nist.gov/cuu/Constants/Table/allascii.txt)。

以下是根据上面的模型公式进行的量纲推导：

\[
\boxed{b_1^{\rm inp}=b_1[\mathrm{fm}]/a,\qquad
b_2^{\rm inp}=b_2[\mathrm{fm}]/a}
\]

\[
\boxed{k_1=\frac{4\pi}{3}\frac{D_1/a^2}{E_{\rm GR}/h}},\qquad
\boxed{k_2=\frac{4\pi}{5}\frac{D_2/a^4}{0.050/h}}.
\]

因此程序输入可直接计算为

\[
k_1=4.07039926151\times10^{-14}\,
\frac{D_1}{E_{\rm GR}[\mathrm{MeV}]},
\]

\[
k_2=1.74427915808\times10^{-22}D_2
=3.48855831616\times10^{-19}S.
\]

当 E1 强度也采用 TRK 近似，可消去 \(D_1\)：

\[
\boxed{k_1=\frac{Z(A-Z)}{2A(m_p/m_e)(E_{\rm GR}/h)^2}}.
\]

另一条等价路线：若已有与本模型归一化一致的 \(\alpha_0^{E1}\)（fm³）或
\(\bar\alpha_0^{E2}\)（fm⁵），分别除以 \(2a^3\)、\(2a^5\)。
核数据表的 \(e^2\mathrm{fm}^{2L}/\mathrm{MeV}\) 尚不是这种体积极化率；
要先乘 \(e^2/(4\pi\epsilon_0)=E_ha_0\)，其 MeV fm 数值为 \(ha\)。
以上两条路线只能选一条，不能再额外乘 \(\alpha_{\rm fs}\)、\(e^2\) 或 \(1/2\)。

## 截断参数的选取与拟合

先定义

\[
A_0=\operatorname{round}(-17.62+2.737Z),\quad
\Delta A=A-A_0,\quad
\delta R_{\rm fm}=R[\mathrm{fm}]-1.2A^{1/3}.
\]

\(R\) 必须和论文的核模型半径同义。不能直接把 rms 电荷半径或者 Fermi
半密度半径代入。选择论文基线 \(R=R_0\) 时明确记录 \(\delta R=0\)；这不是
声称实验半径恰好等于经验半径。

以下为 Eqs. (51), (56) 的 **s 列逐项代入**，输出单位均为 fm：

\[
\begin{aligned}
b_1={}&e^{4.88-0.0123Z+2.89\times10^{-5}Z^2}Z^{-0.178}\\
&+e^{2.05-0.00969Z+2.02\times10^{-4}Z^2}Z^{-1.07}\Delta A\\
&+e^{1.94-0.0190Z+5.90\times10^{-5}Z^2}Z^{-0.198}\delta R_{\rm fm},
\end{aligned}
\]

\[
\begin{aligned}
b_2={}&135-0.600Z\\
&+e^{-0.172-0.0351Z+9.47\times10^{-5}Z^2}\Delta A\\
&+e^{1.93-0.0236Z+5.06\times10^{-5}Z^2}Z^{0.386}\delta R_{\rm fm}.
\end{aligned}
\]

[原文系数：arXiv Tables V–VI](https://arxiv.org/html/2102.00355v1#S4.SS3)。
有对应核素时，也可直接采用原始数值表的 \(b_{ns}\)、\(\tilde b_{ns}\)，记录表号与行。

### 已发现但尚未由作者澄清的文献内部不一致

- Eq. (56) 末项指数印为无波浪号的 \(\tau_3\)，其系数表却给
  \(\tilde\tau_3=0.386\)。上面使用同一 E2 表的 0.386；这是依据表结构和
  数值表核查作出的符号解释。
- E1 s 列 \(\nu_2\) 在预印本和作者上传正式版均印为 \(2.02\times10^{-4}\)。
  逐字代入后，部分高 Z 同位素斜率不复现其数值表，例如 Z=82 的
  \(\partial b_1/\partial A\approx0.1223\) fm，而数值表 A=206、208 均写 26.7 fm。
  不能凭怀疑把系数改成另一个指数；核素参数选择应将此差异留作待核查项。

本轮检索没有找到可核实的正式勘误。这不证明勘误不存在。对精密同位素差分，
应优先复核匹配方程或取得作者确认，不能把拟合误差直接视作物理不确定度。

## 换算例：用论文参数重建 164Dy 模型

这只是量纲与模型的演示，不是已验证的实验输入或生产参数推荐。采用 Z=66、
A=164、论文数值表中的 \(\beta_2=0.348\)、\(b_{ns}=31.6\) fm、
\(\tilde b_{ns}=94.7\) fm。两个截断半径直接取表，不使用上面的近似拟合。
[依据：arXiv Tables IX–X](https://arxiv.org/html/2102.00355v1#S4.SS4)。

按 E1 TRK、E2 集体形变公式、\(R=1.2A^{1/3}\) 以及固定 50 keV 归一化，
代数换算得到（额外数字仅便于复算，不代表物理精度）：

| 中间量 | 数值 |
| --- | --- |
| \(R_0\) | 6.56844441 fm |
| \(E_{\rm GR}\) | 14.1849645 MeV |
| \(B_\downarrow(E2)\) | 11193.0708 \(e^2\mathrm{fm}^4\) |
| 对应 \(B_\uparrow(E2)\) | 5.59653540 \(e^2\mathrm{barn}^2\) |

对应输入片段为：

```fortran
  isotope_z = 66
  isotope_a = 164
  k1 = 3.95214d-14
  b1 = 5.97153d-4
  k2 = 1.95238d-18
  b2 = 1.78957d-3
  parameter_source = 'Flambaum2021 arXiv-v1 Tables IX-X; 164Dy; ns cutoffs; E1 TRK; E2 Bdown from beta2=0.348; Ebar=50keV; R=1.2A^(1/3)fm; illustrative'
```

还需按现行接口提供开关、单位和 `model_id` 等字段。此例没有运行原子计算，
也没有从外部评价库重新取得该核素的实验强度和不确定度。

## 适用性与记录要求

论文主要覆盖中重核；常规区域 \(20\le Z\le98\)，另有超重元素外推。
E2 转动估计主要针对偶偶形变核，典型 \(0.2\lesssim\beta_2\lesssim0.35\)。
本文不把该范围以外的核素或奇 A 核自动视为已经验证。
文献拟合精度描述的是对其参考计算的复现，不是多电子 MCDHF、同位素位移或 HFS
的实验误差保证。[论文适用范围](https://arxiv.org/html/2102.00355v1#S3)。

每组参数至少保存：Z/A、核跃迁方向、强度原值/单位/不确定度/来源、
E1 能量选择、E2 的 `Ebar=50 keV` 约定、半径定义、截断参数的表格或拟合来源、
单位转换常数和最终四个输入值。`parameter_source` 应是这份记录的简短标识。
标量 E2 核极化项不能与静态核四极矩超精细作用或张量极化率混为一项。
