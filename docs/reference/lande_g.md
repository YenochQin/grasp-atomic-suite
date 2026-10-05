# Landé g_J 的理论背景与近似对照

> 文档状态：理论参考草稿、非实施说明。
> 核查日期：2026-10-05；源码基准：`8b298a6`。
> 本文用于概念和 LS 近似对照；约化矩阵元约定及算符系数尚未完成独立文献核查，不作为本仓库数值公式的验收依据。代码流程以[RHFS 报告](../implemented/RHFS_gJ_report.md)为准，程序用法见[仓库 README](../../README_ZH.md)。
> 全部文档状态见[分类索引](../README.md)。

MCDHF/RCI 得到 ASF 后，可通过磁矩算符矩阵元计算朗德因子。下文介绍定义、态混合与 LS 近似，用于理解计算结果。

---

# 1. 基本定义

原子态 \(|\Gamma P J M_J\rangle\) 在弱磁场中的一阶 Zeeman 能级移动可写成

\[
\Delta E = \mu_B\, g_J\, M_J\, B ,
\]

其中

- \(\mu_B\) 是玻尔磁子，
- \(J\) 是总角动量，
- \(M_J\) 是其磁量子数，
- \(g_J\) 就是该原子态的朗德因子。

在相对论原子结构理论中，\(g_J\) 不是简单由 \(L,S,J\) 代入非相对论 Landé 公式得到，而应由 ASF 的波函数直接求出。

---

# 2. MCDHF 中 ASF 的形式

MCDHF 计算得到的原子态函数（ASF）写成 CSF 的线性组合：

\[
|\Gamma P J M_J\rangle
=
\sum_r c_r \, |\gamma_r P J M_J\rangle,
\]

其中

- \(c_r\) 是展开系数，
- \(|\gamma_r P J M_J\rangle\) 是具有相同宇称 \(P\)、总角动量 \(J\) 的组态态函数 CSF。

因此，任何物理量，包括 \(g_J\)，都可以通过该 ASF 对相应算符的矩阵元求得。

---

# 3. 相对论框架下 \(g_J\) 的计算思路

## 3.1 从磁矩定义出发

电子体系的磁矩算符通常写为

\[
\boldsymbol{\mu}
=
-\mu_B \sum_i \left( \mathbf{l}_i + g_s \mathbf{s}_i \right),
\]

其中 \(g_s \approx 2.00232\)；若先忽略 QED 修正，常取 \(g_s=2\)。

对于定态 \(|\Gamma J M_J\rangle\)，朗德因子定义为

\[
g_J
=
-\frac{1}{\mu_B}
\frac{\langle \Gamma J M_J | \mu_z | \Gamma J M_J\rangle}{M_J},
\qquad (M_J\neq 0).
\]

利用 Wigner-Eckart 定理，也可写成与约化矩阵元相关的形式，这是程序实现时更常见的方式。

---

## 3.2 相对论原子结构中的常用公式

在 MCDHF/RCI 框架中，\(g_J\) 常通过磁矩张量算符的一阶矩阵元表示。常见写法为

\[
g_J
=
\frac{2}{\sqrt{J(J+1)(2J+1)}}
\left\langle \Gamma J \left\| \sum_i \left[ -i\frac{\sqrt{2}}{2\alpha^2} r_i
\left( \boldsymbol{\alpha}_i \mathbf{C}^{(1)}_i \right)^{(1)}
+ \frac{g_s-2}{2}\,\beta_i \mathbf{\Sigma}_i
\right] \right\| \Gamma J \right\rangle ,
\]

不过不同文献和程序中算符写法略有不同。对于实际使用者，更重要的是理解：

- **本质上是对总磁偶极矩算符的期望值求解**；
- 在 GRASP/MCDHF 程序中，通常已经内置了 \(g_J\) 的计算模块，不需要你手工写出狄拉克矩阵形式。

---

# 4. 在 LS 耦合近似下的对照公式

如果该态近似是纯的 \(LSJ\) 态，即

\[
|\Gamma J\rangle \approx |\gamma LSJ\rangle,
\]

则可以用熟悉的非相对论 Landé 公式近似：

\[
g_J^{(LS)}
=
1+\frac{J(J+1)+S(S+1)-L(L+1)}{2J(J+1)}.
\]

这在态混合不强时可作为校验。

例如：

- \(^3F_4\): \(L=3, S=1, J=4\)

\[
g_J = 1+\frac{4(5)+1(2)-3(4)}{2\cdot 4\cdot 5}
=1+\frac{20+2-12}{40}
=1.25
\]

但对于像 Ni I 这种**开壳层、强组态混合、显著相对论修正**的体系，单纯用 LS 公式通常不够准确，必须使用 ASF 的实际展开系数来算。

---

# 5. 若 ASF 是多个 LSJ 分量混合，如何由 ASF 求 \(g_J\)

如果 ASF 可以投影到若干个 \(LSJ\) 分量：

\[
|\Gamma J\rangle = \sum_k a_k |\gamma_k L_k S_k J\rangle,
\]

那么一个常见近似是

\[
g_J \approx \sum_k |a_k|^2 g_J(L_kS_kJ),
\]

这里 \(g_J(L_kS_kJ)\) 用普通 Landé 公式给出。

但要注意：

- 这是**近似公式**；
- 严格结果应来自**磁矩算符在 ASF 上的完整矩阵元**；
- 当不同组态、不同 \(LS\) 分量间的非对角贡献不可忽略时，这个加权平均会有偏差。

因此，**正确做法是直接用 MCDHF/RCI 波函数计算磁矩矩阵元**。

---

# 6. 用 ASF 计算 \(g_J\) 的具体步骤

---

## 第一步：完成 MCDHF/RCI 计算，得到收敛 ASF

先得到目标态的

- 轨道函数；
- CSF 展开系数 \(c_r\)；
- 若有必要，包括 Breit 和 QED 修正后的最终波函数。

对于 Ni I 这类复杂体系，建议 \(g_J\) 基于最终的 RCI 波函数而不是仅 MCDHF 波函数，因为：

- 电子关联会改变态混合；
- 态混合对 \(g_J\) 非常敏感。

---

## 第二步：构造磁矩算符矩阵元

在 CSF 基上求

\[
\langle \gamma_r PJ || \mathbf{N}^{(1)} || \gamma_s PJ \rangle
\]

或等价的磁偶极矩算符约化矩阵元，然后由 ASF 展开系数做求和：

\[
\langle \Gamma PJ || \mathbf{N}^{(1)} || \Gamma PJ \rangle
=
\sum_{r,s} c_r c_s
\langle \gamma_r PJ || \mathbf{N}^{(1)} || \gamma_s PJ \rangle.
\]

最后代入 \(g_J\) 的表达式。

---

## 第三步：由约化矩阵元得到 \(g_J\)

常用形式是

\[
g_J
=
\frac{2}{\sqrt{J(J+1)(2J+1)}}
\langle \Gamma J || \mathbf{N}^{(1)} || \Gamma J \rangle.
\]

这里 \(\mathbf{N}^{(1)}\) 表示对应的磁偶极一阶张量算符。具体常数因程序定义可能略有差异，但程序会按其内部约定输出最终 \(g_J\)。

---

# 7. 在 GRASP 中如何得到 \(g_J\)

如果你用的是 GRASP2018/GRASP2K 一类程序，通常流程是：

1. 用 **MCDHF** 优化轨道；
2. 用 **RCI** 加入更充分的关联、Breit、QED；
3. 用相应的属性模块计算 \(g_J\)、超精细常数、跃迁矩阵元等。

不同版本模块名略有区别，但本质上是：

- 读取 `.c`/`.w` 等波函数文件；
- 指定要计算的态；
- 调用 property 程序求磁矩相关量；
- 输出每个 ASF 的 \(g_J\)。

本仓库的具体程序名、输入文件和执行方式见[仓库 README](../../README_ZH.md)。

---

# 8. 如何判断计算得到的 \(g_J\) 是否可靠

对 Ni I 这类复杂体系，\(g_J\) 是一个非常好的**态标识诊断量**。可从以下几方面判断：

## 8.1 与实验值比较

若有实验 Zeeman 数据，比较最直接。
\(g_J\) 往往比能量对态混合更敏感，因此：

- 能量接近并不代表态标识对；
- \(g_J\) 对了，往往说明波函数成分更可信。

---

## 8.2 与纯 LS 值比较

若某态理论 \(g_J\) 明显偏离 Landé 公式值，可能意味着：

- 态存在显著中间耦合；
- 组态混合很强；
- 或态排序有问题。

---

## 8.3 检查态混合与能级交叉

Ni I 中容易出现近简并态。若两个同 \(J^\pi\) 态靠得很近，则：

- 微小的关联模型变化
- Breit 修正
- 轨道优化策略

都可能改变混合系数，从而显著改变 \(g_J\)。

所以做 \(g_J\) 时，建议同时检查：

- ASF 主导组态百分比；
- \(LSJ\) 组成；
- 与相邻态的相互作用强弱。

---

# 9. 一个简单实用的理解公式

如果你暂时不考虑严格相对论算符，而只想从 MCDHF 输出的态成分估算 \(g_J\)，可以这样做：

1. 先用 `jj2lsj` 或类似分析工具，把 ASF 分解到主要 \(LSJ\) 成分；
2. 对每个主要成分算 Landé 因子：

\[
g_k = 1+\frac{J(J+1)+S_k(S_k+1)-L_k(L_k+1)}{2J(J+1)};
\]

3. 用权重 \(w_k\) 加权：

\[
g_J^{\text{approx}} = \sum_k w_k g_k,
\qquad \sum_k w_k \approx 1.
\]

这可以帮助你快速理解某个 ASF 的 \(g_J\) 来源，但不要替代正式结果。

---

# 10. 结论

用 MCDHF 方法得到 ASF 后，朗德因子 \(g_J\) 的计算核心是：

\[
\boxed{
g_J
=
-\frac{1}{\mu_B}
\frac{\langle \Gamma J M_J | \mu_z | \Gamma J M_J\rangle}{M_J}
}
\]

或者等价地，通过磁偶极矩张量算符的约化矩阵元计算。
在 MCDHF/RCI 中，具体就是：

- 用 ASF 的 CSF 展开系数；
- 计算磁矩算符在 CSF 基底上的矩阵元；
- 再组合得到态的 \(g_J\)。
