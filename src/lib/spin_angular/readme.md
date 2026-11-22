# 自旋–角向系数

`spin_angular` 计算两个 $jj$ 耦合组态态函数（CSF）之间的单电子和双电子自旋–角向系数。以下公式对应 suite `2d4c721` 加 `spin_angular` 命名迁移后的实现。文件、模块与旧例程名的对应关系见[命名对照](../../../docs/implemented/spin_angular_naming.md)，形参与工作数组见[变量对照](../../../docs/implemented/librang90_variable_migration.md)。

## 入口和符号

| 文件 | 例程 | 输出 |
| --- | --- | --- |
| [one_electron_scalar_coefficients.f90](one_electron_scalar_coefficients.f90) | `calculate_one_electron_scalar_coefficients` | 单电子偶标量算符的角向系数 |
| [one_electron_tensor_coefficients.f90](one_electron_tensor_coefficients.f90) | `calculate_one_electron_tensor_coefficients` | 给定秩 $K$ 和宇称的单电子张量系数 |
| [two_electron_angular_coefficients.f90](two_electron_angular_coefficients.f90) | `calculate_two_electron_angular_coefficients` | 按 CSF 对、四个轨道标签和阶数输出 Coulomb 或 Breit 角向系数 |

$L$ 和 $R$ 分别标记 bra 和 ket。$j$ 是单电子角动量，$J$ 是子壳层项、耦合链节点或整个 CSF 的总角动量，$Q$ 是准自旋。本文使用 $[x]=2x+1$；源码用 `two_j_1` 等整数存储 $2j$，`two_j_plus_one` 存储 $[j]$。活跃子壳层的位置通过 `peel_orbital_indices` 转换为完整轨道表中的编号。

为对应源码的整除测试，定义 $\eta_p(u)=1$ 当 $u\equiv0\pmod p$，其余情况 $\eta_p(u)=-1$。$u$ 是源码里的整数相位量。例如 $u$ 为偶数时，$\eta_4(u)=(-1)^{u/2}$。

令 $\Phi_L,\Phi_R$ 为两个 CSF。以调用方使用的径向积分为基准，Dirac–Coulomb 矩阵元的装配形式为

$$
H_{LR}^{\mathrm{DC}}=
\sum_{ab}T_{LR}(ab)\,I(ab)
+\sum_{abcd,k}V_{LR}^{k}(abcd)\,R^k(abcd).
$$

$T$ 对应 `shell_coefficients` 输出的单电子权重，$V^k$ 对应传给 `speak` 的 Coulomb 角向贡献。$I$ 是单电子积分，$R^k$ 是 Slater 双电子径向积分。Coulomb 的 $V^k$ 已含单轨道球张量因子。Breit 路径传给 `collect_breit_angular_coefficient` 的权重按径向积分类型累加。见 [RMCDHF 的矩阵装配](../../apps/rmcdhf/assemble_hamiltonian_block.f90)和 [RCI 的 Coulomb/Breit 装配](../../apps/rci/assemble_hamiltonian_block.f90)。

## 三角选择、Racah 与 Wigner 符号

### 三角条件和共同阶数区间

[angular_triangle_and_parity.f90](angular_triangle_and_parity.f90) 检查整数标签 $a,b,c$ 是否满足 $|a-b|\le c\le a+b$ 且 $a+b+c$ 为偶数。用于角动量耦合时传入二倍角动量；Coulomb 路径用整数 $l,l',k$，同时实施球谐张量的宇称选择。

[shared_tensor_rank_range.f90](shared_tensor_rank_range.f90) 对两对二倍角动量 $(a,b)$、$(c,d)$ 求公共区间：

$$
\lambda_{\min}=\max(|a-b|,|c-d|),\qquad
\lambda_{\max}=\min(a+b,c+d).
$$

返回 $\lambda_{\min}$ 和 $\lambda_{\max}-\lambda_{\min}+1$。双电子调用方在此基础上按步长 2 遍历允许的阶数；循环变量常采用 $2k+1$。

### Racah 有限和与 $6j$

对符号 $\begin{Bmatrix}a&b&c\\d&e&f\end{Bmatrix}$，
[wigner_6j_triangle_selection.f90](wigner_6j_triangle_selection.f90) 检查
$(a,b,c)$、$(a,e,f)$、$(d,b,f)$ 和 $(d,e,c)$ 四组三角。定义

$$
\begin{aligned}
\Delta(x,y,z)&=\left[\frac{(x+y-z)!(x-y+z)!(-x+y+z)!}{(x+y+z+1)!}\right]^{1/2},\\
\alpha&=(a+b+c,\ a+e+f,\ d+b+f,\ d+e+c),\\
\beta&=(a+b+d+e,\ a+c+d+f,\ b+c+e+f).
\end{aligned}
$$

用 $z_0=\max_i\alpha_i$ 和 $z_1=\min_i\beta_i$，Racah 有限和为

$$
\begin{Bmatrix}a&b&c\\ d&e&f\end{Bmatrix}
=\Delta(a,b,c)\Delta(a,e,f)\Delta(d,b,f)\Delta(d,e,c)
\sum_{z=z_0}^{z_1}\frac{(-1)^z(z+1)!}
{\prod_{i=1}^{4}(z-\alpha_i)!\prod_{i=1}^{3}(\beta_i-z)!}.
$$

[racah_coefficients.f90](racah_coefficients.f90) 和 [racah_coefficients_scott_fischer.f90](racah_coefficients_scott_fischer.f90) 使用上式的相邻项比值：

$$
\rho_z=
\frac{(z+1)\prod_{i=1}^{3}(\beta_i-z+1)}
{\prod_{i=1}^{4}(z-\alpha_i)},\qquad
h_{z_1}=1,\qquad h_{z-1}=1-\rho_z h_z.
$$

实现从上界递推至下界，再乘上带符号的首项与四个三角因子。阶乘由 $\log(n!)$ 表计算，表的第 $n+1$ 项存储 $\log(n!)$。两例程的参数按 $W(a,b,e,d;c,f)$ 排列。与 Wigner 符号的关系为

$$
W(a,b,e,d;c,f)=(-1)^{a+b+e+d}
\begin{Bmatrix}a&b&c\\d&e&f\end{Bmatrix}.
$$

[wigner_6j.f90](wigner_6j.f90) 在通用路径乘回这个相位。六个二倍参数均在 $0$ 到 $4$ 时使用小表缓存。零角动量路径直接用相位和简并度。以 $a=0$ 为例，允许的三角强制 $b=c$、$e=f$，此时值为 $(-1)^{d+e+c}/\sqrt{[e][c]}$。标准有限和参见 [DLMF `34.4](https://dlmf.nist.gov/34.4#E2)。

### 小角动量的 $6j$ 路径

当最小角动量为 $\tfrac12,1,\tfrac32,2,3,4$ 时，`sixj` 通过列交换和成对上下交换将该值移到第六格 $f=t$，再使用以下实现。$t$ 是角动量值，`two_j_*` 中存的量为 $2t$。

| $t$ | 源文件 | 计算 |
| --- | --- | --- |
| $\tfrac12$ | [wigner_6j_with_j_one_half.f90](wigner_6j_with_j_one_half.f90) | 按两条边差的符号选根式 |
| $1$ | [wigner_6j_with_j_one.f90](wigner_6j_with_j_one.f90) | 按两条边差的符号和相等关系选闭式 |
| $\tfrac32$ | [wigner_6j_with_j_three_halves.f90](wigner_6j_with_j_three_halves.f90) | 用二倍边差 $\pm1,\pm3$ 选闭式 |
| $2$ | [wigner_6j_with_j_two.f90](wigner_6j_with_j_two.f90) | 用二倍边差 $0,\pm2,\pm4$ 选闭式 |
| $3$ | [wigner_6j_with_j_three.f90](wigner_6j_with_j_three.f90) | 用二倍边差 $0,\pm2,\pm4,\pm6$ 选闭式 |
| $4$ | [wigner_6j_with_j_four.f90](wigner_6j_with_j_four.f90) | 用 $t=2,3$ 的值递推，或回退至 Racah 和 |

前五行的闭式可以用同一个短有限和来理解。将 $f=t$ 代入 $\alpha,\beta,\Delta$，把首项提到求和号外，则乘在首项上的因子为

$$
\sum_{r=0}^{z_1-z_0}(-1)^r\prod_{u=1}^{r}\rho_{z_0+u}.
$$

空积为 1。三角条件使 $z_1-z_0\le 2t$，所以各小 $t$ 分支只包含有限的多项式与根式。源码以 $2a-2e$ 和 $2b-2d$ 的具体差值选择这些分支。以 $t=\tfrac12$、$b<d$、$a<e$ 为例，令 $s=a+b+c$，`sixj5` 的计算为

$$
S_{1/2}=(-1)^{s+1}
\sqrt{\frac{(s+2)(s-2c+1)}
{(2b+1)(2b+2)(2a+1)(2a+2)}}.
$$

对于 $t=1,\tfrac32,2,3$，若 $a-e=b-d=t$，这四个例程的相应分支可用下降阶乘统一表示。令 $n=2t$，$(x)_{\downarrow m}=\prod_{v=0}^{m-1}(x-v)$，则

$$
S_t=(-1)^s
\sqrt{
\frac{(s+1)_{\downarrow n}\,(s-2c)_{\downarrow n}}
     {(2a+1)_{\downarrow(n+1)}(2b+1)_{\downarrow(n+1)}}
}.
$$

这里四个例程分别取 $n=2,3,4,6$。在 $t=1$ 的 $a=e,\ b=d$ 分支，令 $A=2c,\ B=2b,\ C=2a$，闭式改为

$$
S_1=(-1)^{s+1}
\frac{(B^2+C^2-A^2)/2+B+C-A}
     {\sqrt{B(B+1)(B+2)C(C+1)(C+2)}}.
$$

对 $t=4$，定义 $S_t=\begin{Bmatrix}a&b&c\\d&e&t\end{Bmatrix}$，$u=3$，$g(x)=x(x+1)$，以及

$$
\begin{aligned}
F_+(x,y)&=(u+x+y+2)(u-x+y+1)(u+x-y+1)(-u+x+y),\\
F_-(x,y)&=(u+x+y+1)(u-x+y)(u+x-y)(-u+x+y+1),\\
X_+&=u\sqrt{F_+(d,b)F_+(a,e)},\\
X_-&=(u+1)\sqrt{F_-(d,b)F_-(a,e)},\\
X_0&=[u]\bigl[2\{g(u)g(c)-g(d)g(a)-g(b)g(e)\}\\
&\hspace{4em}+\{g(u)-g(d)-g(b)\}\{g(u)-g(a)-g(e)\}\bigr].
\end{aligned}
$$

若 $t=3$ 的三角条件成立，`sixj4` 取 $S_4=(X_0S_3-X_-S_2)/X_+$。两个分子权重的绝对值小于代码容差时直接跳过对应项。$t=3$ 的三角条件不成立时由 Racah 有限和计算 $S_4$。

### $9j$ 的乘积和与零参数

[wigner_9j.f90](wigner_9j.f90) 的筛选模式检查三行和三列的三角关系，求值模式使用

$$
\begin{Bmatrix}a&b&c\\d&e&f\\g&h&i\end{Bmatrix}=
\sum_x(-1)^{2x}[x]
\begin{Bmatrix}a&b&c\\f&i&x\end{Bmatrix}
\begin{Bmatrix}d&e&f\\b&x&h\end{Bmatrix}
\begin{Bmatrix}g&h&i\\x&a&d\end{Bmatrix},
$$

其中 $\max(|a-i|,|f-b|,|d-h|)\le x\le\min(a+i,f+b,d+h)$。实现按步长 2 遍历 $2x$。例如 $a=0$ 时，[wigner_9j_with_zero.f90](wigner_9j_with_zero.f90) 在相容的三角条件下返回

$$
\begin{Bmatrix}0&b&c\\d&e&f\\g&h&i\end{Bmatrix}=
\frac{(-1)^{h+b+f+d}}{\sqrt{[b][d]}}
\begin{Bmatrix}e&h&b\\i&f&d\end{Bmatrix}.
$$

其余零位置使用同样的一次 $6j$ 求值与对应的相位和简并度。参见 [DLMF `34.6](https://dlmf.nist.gov/34.6#E2) 和 [`34.7](https://dlmf.nist.gov/34.7)。

## 子壳层项、准自旋和张量矩阵元

### 工作数组与项索引

[jj_subshell_state_packing.f90](jj_subshell_state_packing.f90) 用 Dirac 角量子数 $\kappa$、占据数 $N$ 和 seniority $\nu$ 构造

$$
2j=2|\kappa|-1,\qquad
l=\frac{2j+\operatorname{sgn}\kappa}{2},\qquad
2Q=|\kappa|-\nu,\qquad
M_Q=\frac{N-|\kappa|}{2}.
$$

每侧整数数组的七个槽依次存储项索引、$n$、$2j$、$N$、$l$、$2J$、$2Q$，实数数组的三个槽存储 $Q,J,M_Q$。[jj_active_subshell_states.f90](jj_active_subshell_states.f90) 为至多四个活跃子壳层填充这些数组。

[jj_subshell_term_lookup.f90](jj_subshell_term_lookup.f90) 用 $2\nu=2j+1-2(2Q)$ 构造基数 100 的键：

$$
B=10^4(2Q)+10^2(2\nu)+2J.
$$

[jj_subshell_term_decoding.f90](jj_subshell_term_decoding.f90) 通过 [packed_integer_digit.f90](packed_integer_digit.f90) 反解。第 $p$ 个字段为 $\lfloor B/100^{p-1}\rfloor\bmod100$，$p=1,2,3$。少电子 $j=9/2$ 表的项索引带有 300 的偏移。

[quasispin_projection_selection.f90](quasispin_projection_selection.f90) 检查 $|M_Q|\le Q$ 且 $2Q+2|M_Q|$ 为偶数。[同奇偶](jj_intermediate_terms_same_parity.f90)、[bra 奇偶](jj_intermediate_terms_bra_parity.f90)、[ket 奇偶](jj_intermediate_terms_ket_parity.f90)三个中间项选择例程用两侧的项表范围和角动量三角关系筛选。若 ket 到中间项的占据变化为整数 $\delta$，传入求和的中间投影为 $M_{Q,I}=M_{Q,R}+\delta/2$。

[subshell_tensor_diagnostics.f90](subshell_tensor_diagnostics.f90) 报告不支持的高 $j$ 占据或项编码，并终止该路径。

### 准自旋 Clebsch–Gordan 系数

以下 $C=\langle Q_R M_R;k_q m\mid Q_L M_L\rangle$，用 $x=Q_L$、$z=M_L$。所有分支先要求 $M_L=M_R+m$、投影在范围内以及 $(Q_R,k_q,Q_L)$ 可构成三角。

[quasispin_cg_rank_half.f90](quasispin_cg_rank_half.f90) 取 $s=2m=\pm1$：

$$
C_{\frac12,m}=
\begin{cases}
\sqrt{(x+sz)/(2x)},& x=Q_R+\tfrac12,\\
-s\sqrt{(x-sz+1)/(2x+2)},& x=Q_R-\tfrac12.
\end{cases}
$$

[quasispin_cg_rank_one_zero_projection.f90](quasispin_cg_rank_one_zero_projection.f90) 使用

$$
C_{1,0}=
\begin{cases}
\sqrt{(x+z)(x-z)/\{x(2x-1)\}},&x=Q_R+1,\\
z/\sqrt{x(x+1)},&x=Q_R,\\
-\sqrt{(x+z+1)(x-z+1)/\{(x+1)(2x+3)\}},&x=Q_R-1.
\end{cases}
$$

当 $M_R=M_L=0$ 时还检查 $Q_R+Q_L+1$ 的奇偶性。[quasispin_cg_rank_one_unit_projection.f90](quasispin_cg_rank_one_unit_projection.f90) 对 $m=s=\pm1$ 使用

$$
C_{1,s}=
\begin{cases}
\sqrt{(x+sz-1)(x+sz)/\{(2x-1)(2x)\}},&x=Q_R+1,\\
-s\sqrt{(x+sz)(x-sz+1)/\{2x(x+1)\}},&x=Q_R,\\
\sqrt{(x-sz+1)(x-sz+2)/\{(2x+2)(2x+3)\}},&x=Q_R-1.
\end{cases}
$$

[quasispin_cg_zero_projection.f90](quasispin_cg_zero_projection.f90) 对 $k_q=0$ 返回 $\delta_{Q_LQ_R}\delta_{M_LM_R}$，秩 1 时转给 $C_{1,0}$ 路径。

### 单张量 $A$ 和双张量 $W$

$A$ 表示一个产生或湮灭张量，$W$ 表示两个 $A$ 耦合成的张量。$\widetilde A$ 和 $\widetilde W$ 同时对角动量与准自旋约化。普通表值由带符号的整数分子 $n$ 和分母 $d$ 还原：

$$
\operatorname{sgn}(n)\sqrt{\frac{|n|}{d}}.
$$

[jj_reduced_cfp.f90](jj_reduced_cfp.f90) 依 $j$ 选择表并处理 bra、ket 对调。低 $j$ 表包括 $j=\tfrac12$ 的常数 $-2$ 和 $j=\tfrac32$ 的负根式项。对照低 $j$ 或 `jj_reduced_cfp_j_nine_halves` 查表路径，若两侧项索引逆序，代码用以下整数相位调整：

$$
\eta_4(2Q_L+2J_L-2Q_R-2J_R+2j-1).
$$

[jj_reduced_cfp_j_nine_halves.f90](jj_reduced_cfp_j_nine_halves.f90) 使用相应的固定或逐项分母。[jj_reduced_cfp_two_electron_terms.f90](jj_reduced_cfp_two_electron_terms.f90) 对少电子项从角向矩阵元反解完全约化值。以 $N_L=2$、$N_R=1$ 排列时，

$$
\widetilde A
=\eta_2(-N_L-1)\,\frac{\sqrt{N_L[Q_L][J_L]}}{C_{\frac12,\frac12}}.
$$

用于小 $j$ 张量乘积的单壳层量是 $A_{LR}(m)=C_{\frac12,m}\widetilde A_{LR}/\sqrt{[Q_L]}$。在高 $j$ 的少电子路径，[jj_a_matrix_element_few_electrons.f90](jj_a_matrix_element_few_electrons.f90) 直接使用

$$
A_{LR}(m)=
\begin{cases}
\eta_4(2J_R+2j-2J_L+2N_R)\sqrt{[J_R]N_R},& m=-\tfrac12,\\
\eta_2(N_L)\sqrt{[J_L]N_L},& m=+\tfrac12 .
\end{cases}
$$

[jj_reduced_w.f90](jj_reduced_w.f90) 按 $(j,k_q,k)$ 选择 $\widetilde W^{(k_q,k)}$ 的表值，$k_q$ 为准自旋秩，$k$ 为角动量秩。表分为 [$j=\tfrac12$](jj_reduced_w_j_one_half.f90)、[$\tfrac32$](jj_reduced_w_j_three_halves.f90)、[$\tfrac52$](jj_reduced_w_j_five_halves.f90)、[$\tfrac72$ 低秩](jj_reduced_w_j_seven_halves.f90)和 [$\tfrac72$ 高秩](jj_reduced_w_j_seven_halves_higher_ranks.f90)。最后一组在调换项索引时乘 $\eta_4(2J_R-2J_L+2Q_R-2Q_L)$。

[jj_reduced_w_from_cfp.f90](jj_reduced_w_from_cfp.f90) 给出同一类完全约化量的求和计算：

$$
\begin{aligned}
\widetilde W_{LR}^{(k_q,k)}={}&
\eta_4(2Q_L+2J_L+2Q_R+2J_R+2k_q+2k)\sqrt{[k_q][k]}\\
&{}\times \sum_I
\begin{Bmatrix}j&j&k\\J_R&J_L&J_I\end{Bmatrix}
\begin{Bmatrix}\tfrac12&\tfrac12&k_q\\Q_R&Q_L&Q_I\end{Bmatrix}
\widetilde A_{LI}\widetilde A_{IR}.
\end{aligned}
$$

[jj_w_matrix_element.f90](jj_w_matrix_element.f90) 用两个 $A$ 的准自旋投影 $m_1,m_2$ 恢复占据依赖。令 $m=m_1+m_2$，通过对应的选择规则后，

$$
W_{LR}^{k}(m_1,m_2)=
\begin{cases}
C_{1,m}\widetilde W_{LR}^{(1,k)}/\sqrt{[Q_L]},&
m\ne0,\quad k\text{ 偶},\quad N_L-N_R=2m,\\
([j]-N_R)\sqrt{[J_L]/[j]},&m=0,\quad k=0,\quad m_1<0,\\
-N_R\sqrt{[J_L]/[j]},&m=0,\quad k=0,\quad m_1>0,\\
\sigma\, C_{k_q,0}\widetilde W_{LR}^{(k_q,k)}/\sqrt{2[Q_L]},&
m=0,\quad k\ne0 .
\end{cases}
$$

$m=0$ 的分支要求 $N_L=N_R$，$k=0$ 还要求项索引相同。最后一行在 $k$ 偶时取 $k_q=1$，奇时取 $k_q=0$；$k$ 为奇且 $m_1<0$ 时 $\sigma=-1$，其余 $\sigma=1$。

[jj_w_matrix_element_high_j.f90](jj_w_matrix_element_high_j.f90) 用少电子 $A$ 直接求和：

$$
W_{LR}^{k}
=\eta_4(2k+2J_L+2J_R)\sqrt{[k]}
\sum_I\begin{Bmatrix}j&j&k\\J_R&J_L&J_I\end{Bmatrix}
A_{LI}(m_1)A_{IR}(m_2).
$$

求和中的中间项占据从 ket 侧加上 $2m_2$ 后获得。各张量约化的理论背景与准自旋约定见 [Gaigalas, Atoms 10, 129, ``2.3 和 3.3](https://epublications.vu.lt/object/elaba%3A144935795/144935795.pdf)。

### 同一子壳层的 $AW$、$WA$ 和 $WW$

令 $K$ 为复合张量的秩，$k$ 为 $W$ 的秩，$P=\eta_4(2K+2J_L+2J_R+\zeta)$。下面两式在小 $j$ 的 [AW](jj_aw_matrix_element.f90) 和 [WA](jj_wa_matrix_element.f90) 路径取 $\zeta=2$，并使用含 Clebsch–Gordan 系数的 $A$。在高 $j$ 的 [AW](jj_aw_matrix_element_high_j.f90) 和 [WA](jj_wa_matrix_element_high_j.f90) 路径取 $\zeta=0$，$A$ 来自少电子公式。

$$
\begin{aligned}
(AW)_{LR}^{K}
&=P\sqrt{[K]}\sum_I
\begin{Bmatrix}j&k&K\\J_R&J_L&J_I\end{Bmatrix}
A_{LI}(m_1)W_{IR}^{k}(m_2,m_3),\\
(WA)_{LR}^{K}
&=P\sqrt{[K]}\sum_I
\begin{Bmatrix}k&j&K\\J_R&J_L&J_I\end{Bmatrix}
W_{LI}^{k}(m_1,m_2)A_{IR}(m_3).
\end{aligned}
$$

$WA$ 的中间占据为 $N_R+2m_3$。小 $j$ 的 $AW$ 按 $N_R+2m_2+2m_3$ 选中间项。高 $j$ 的 $AW$ 源码在 `intermediate_occupation_change` 中使用的是 $N_R+4m_2$。

[jj_ww_scalar_matrix_element.f90](jj_ww_scalar_matrix_element.f90) 检查 $J_L=J_R$，将两个同秩 $W$ 耦合为标量：

$$
(WW)_{LR}^{0}
=\frac{1}{\sqrt{[k][J_L]}}\sum_I
\eta_4(2k-2J_L+2J_I)
W_{LI}^{k}(m_1,m_2)W_{IR}^{k}(m_3,m_4).
$$

中间项由 $N_R+2m_3+2m_4$ 和相应准自旋投影确定。

### 多个子壳层的张量乘积

下表中的 $A_s$、$W_s$ 等是第 $s$ 个子壳层上已包含准自旋投影依赖的角向矩阵元。$A$ 改变占据 $2m$，$W$ 改变占据 $2(m_1+m_2)$，$AW$ 或 $WA$ 改变占据 $2(m_1+m_2+m_3)$。各文件先检查实际占据差，再计算乘积。

| 文件 | 返回的子壳层乘积 $M$ |
| --- | --- |
| [subshell_product_a_a.f90](subshell_product_a_a.f90) | $A_1A_2$ |
| [subshell_product_w_w.f90](subshell_product_w_w.f90) | $W_1W_2$ |
| [subshell_product_wa_a.f90](subshell_product_wa_a.f90) | $-(WA)_1 A_2$ |
| [subshell_product_a_aw.f90](subshell_product_a_aw.f90) | $-A_1(AW)_2$ |
| [subshell_product_a_a_w.f90](subshell_product_a_a_w.f90) | $A_1A_2W_3$ |
| [subshell_product_a_a_a_a.f90](subshell_product_a_a_a_a.f90) | $A_1A_2A_3A_4$ |

## 耦合链重耦合

### 单个节点

以 $x_L,x_R$ 表示节点处受张量影响的输入角动量，$s$ 表示其余耦合伙伴的角动量，$y_L,y_R$ 表示耦合后的角动量。[first_active_shell_recoupling.f90](first_active_shell_recoupling.f90) 的节点因子为

$$
D_1(k)=
\eta_4(2s+2y_L+2x_R+2k)
\sqrt{[x_L][y_R]}
\begin{Bmatrix}k&x_R&x_L\\s&y_L&y_R\end{Bmatrix}.
$$

活跃位置为第 1 壳层时，再乘交换首对耦合顺序的相位 $\eta_4(2x_L+2x_R+4s-2y_L-2y_R)$。[final_shell_recoupling.f90](final_shell_recoupling.f90) 的 $D_5$ 也使用这个节点表达式。活跃位置在末壳层时，输入是末壳层的两侧角动量，$s$ 是前缀耦合；活跃位置在末壳层前面时，输入是两侧前缀耦合，$s$ 是末壳层角动量，并乘上述交换相位。

[spectator_shell_chain_recoupling.f90](spectator_shell_chain_recoupling.f90) 返回中间旁观节点的乘积 $D_3(i,j;k)$。对于每个这样的节点，$x_L,x_R$ 是前缀耦合，节点因子为

$$
\eta_4(2k+2s+2x_L+2y_R)
\sqrt{[x_L][y_R]}
\begin{Bmatrix}k&x_R&x_L\\s&y_L&y_R\end{Bmatrix}.
$$

乘积到位置 $j$ 之前结束。从位置 $i=1$ 开始时跳过首对壳层共享的耦合节点。

[later_active_shell_recoupling.f90](later_active_shell_recoupling.f90) 和 [pair_tensor_recoupling.f90](pair_tensor_recoupling.f90) 使用另一种配置。记 $x_L,x_R$ 为前缀耦合，$s_L,s_R$ 为后加入的壳层角动量，$y_L,y_R$ 为耦合输出：

$$
\begin{aligned}
D_2(k)&=\eta_4(2y_L+2x_R+2s_L+2k)
   \sqrt{[s_R][x_L]}
   \begin{Bmatrix}k&s_R&s_L\\y_L&x_L&x_R\end{Bmatrix},\\
D_4(k_1,k_2;K)&=\sqrt{[x_L][K][s_L][y_R]}
   \begin{Bmatrix}x_R&k_1&x_L\\s_R&k_2&s_L\\y_R&K&y_L\end{Bmatrix}.
\end{aligned}
$$

`evaluate_coefficient=0` 时，节点例程用同一组 $6j$ 或 $9j$ 的三角条件筛选。候选 CSF 的旁观壳层量子数与中间 $J$ 另由[标量相互作用](scalar_interaction_recoupling_selection.f90)、[单电子标量](one_electron_scalar_recoupling_selection.f90)和[单电子张量](one_electron_tensor_recoupling_selection.f90)筛选例程检查。

### 两到四个活跃位置

令活跃位置满足 $i_1<i_2<\cdots$，$Z_n=1/\sqrt{\prod_{s=1}^{n}[J_{L,i_s}]}$。为简化表达，$D_r(i,j;k)$ 的位置参数指对应节点例程读取的壳层。在两处最早活跃位置定义

$$
F_{12}(k)=D_1(i_1;k)D_3(i_1,i_2;k).
$$

当 $(i_1,i_2)=(1,2)$ 时取 $F_{12}=1$。当两者之间没有需要传播的节点时，$D_3=1$。其余重耦合因子依次相乘：

- [two_subshell_recoupling.f90](two_subshell_recoupling.f90) 的秩 0 求值直接返回 $Z_2$。秩 $K\ne0$ 的结果为

  $$
  R_2=
  Z_2\sqrt{\frac{[J_{L,i_2}]}{[K][J_{R,i_2}]}}
  D_2(i_1,i_2;K)F_{12}(K).
  $$

- [three_subshell_recoupling.f90](three_subshell_recoupling.f90) 将前两个张量 $(k_1,k_2)$ 耦合至 $K$，再连接第三个位置：

  $$
  \begin{aligned}
  R_3={}&Z_3\sqrt{\frac{[J_{L,i_3}]}{[K][J_{R,i_3}]}}\,
  D_3(i_2,i_3;K)D_2(i_1,i_3;K)\\
  &{}\times D_4(i_1,i_2;k_1,k_2;K)F_{12}(k_1).
  \end{aligned}
  $$

- [four_subshell_recoupling.f90](four_subshell_recoupling.f90) 传递第一配对的秩 $k_{12}$，再将它与第三个张量耦合到 $k_{34}$：

  $$
  \begin{aligned}
  R_4={}&Z_4\sqrt{\frac{[J_{L,i_4}]}{[k_{34}][J_{R,i_4}]}}\,
  D_3(i_2,i_3;k_{12})D_3(i_3,i_4;k_{34})D_2(i_1,i_4;k_{34})\\
  &{}\times D_4(i_1,i_2;k_1,k_2;k_{12})
  D_4(i_2,i_3;k_{12},k_3;k_{34})F_{12}(k_1).
  \end{aligned}
  $$

[three_subshell_ordering.f90](three_subshell_ordering.f90) 返回三个互异位置的最小、中间和最大值。
[permuted_three_subshell_recoupling.f90](permuted_three_subshell_recoupling.f90) 将输入位置排序后调用 $R_3$，并按下表乘相位。数字表示排序后的原参数位置，$r_s=2k_s$。

| 排序 | 相位 | 排序 | 相位 |
| --- | --- | --- | --- |
| 123 | $1$ | 213 | $\eta_4(r_1+r_2-r_3)$ |
| 312 | $\eta_2(r_3)$ | 321 | $\eta_4(r_1+r_2+r_3)$ |
| 132 | $\eta_4(r_1-r_2-r_3)$ | 231 | $\eta_2(r_1)$ |

### 单电子重耦合与输出系数

[one_electron_same_subshell_recoupling.f90](one_electron_same_subshell_recoupling.f90)
从 $Z_1$ 出发，作用位置为 $i$、末端为 $n$ 时按分支乘 $D_5(n,i;K)$、$D_1(i;K)$ 和 $D_3(i,n;K)$。$K=0$ 的求值直接使用 $Z_1$。$n=2$ 时由 $D_1$ 处理首对，$n>2$ 且 $i=n$ 时只需末端节点，$n=1$ 时无需节点因子。

[one_electron_two_subshell_recoupling.f90](one_electron_two_subshell_recoupling.f90)
的非平凡分支为

$$
R^{(1e)}_2=Z_2\,D_3(i_2,n;K)D_5(n,i_2;K)D_4(i_1,i_2;j_{i_1},j_{i_2};K)F_{12}(j_{i_1}).
$$

若 $i_2=n$，末端因子省略。这里 $j_{i_1}$ 和 $j_{i_2}$ 是两个作用轨道的单电子角动量。

[单电子标量入口](one_electron_scalar_coefficients.f90)和[一般张量入口](one_electron_tensor_coefficients.f90)
用 $\sum_s |N_{L,s}-N_{R,s}|$ 识别同占据或单次转移。一般张量要求 $(J_L,K,J_R)$ 能组成三角，并按 `operator_parity` 检查宇称。占据相同的路径把各系数按实际轨道编号写入 `shell_coefficients`；单次转移路径把系数写在第 1 槽，轨道对单独返回。

对同一个壳层，标量入口的该分支只接受同一 CSF。[标量同壳公式](one_electron_scalar_same_subshell.f90)和[一般张量同壳公式](one_electron_tensor_same_subshell.f90)分别为

$$
T^{(0)}_{aa}=-W^0(+\tfrac12,-\tfrac12)\sqrt{\frac{[j_a]}{[J_{L,a}]}},
\qquad
T^{(K)}_{aa}=-W^K(+\tfrac12,-\tfrac12)R^{(1e)}_1\sqrt{\frac{[j_a]}{[K]}}.
$$

对从轨道 $b$ 转至 $a$ 的项，设 $F=A_a(+\tfrac12)A_b(-\tfrac12)$。
$N_{\mathrm{span}}$ 是按活跃壳层序列从较早位置到较晚位置之前的 bra 占据数之和，定义 $\sigma_N=\eta_2(N_{\mathrm{span}})$。
两套分支分别计算：

- [one_electron_scalar_between_subshells.f90](one_electron_scalar_between_subshells.f90) 要求 $j_a=j_b$。令 $\sigma_l=\eta_4(\sum_{s=a,b}l_s(N_{L,s}-N_{R,s}))$，则

  $$
  T^{(0)}_{ab}=\sigma_N\sigma_l F R_2(j_b)\sqrt{[j_a]}.
  $$

- [one_electron_tensor_between_subshells.f90](one_electron_tensor_between_subshells.f90) 的结果为

  $$
  T^{(K)}_{ab}=\sigma_N \sigma_{\mathrm{ord}} F R^{(1e)}_2\sqrt{\frac{[j_a]}{[K]}},
  \qquad
  \sigma_{\mathrm{ord}}=
  \begin{cases}
  \eta_4(2j_a+2j_b-2K+2),&a\text{ 晚于 }b,\\
  1,&a\text{ 早于 }b.
  \end{cases}
  $$

## 双电子系数

### 按占据差分派

[two_electron_angular_coefficients.f90](two_electron_angular_coefficients.f90) 要求两侧 CSF 总 $J$ 和宇称相同。设 $d_s=N_{L,s}-N_{R,s}$，使用

$$
D=\sum_s|d_s|,\qquad
G=\#\{s:|d_s|=1\}+20\,\#\{s:|d_s|=2\}.
$$

| 条件 | 双电子角向路径 |
| --- | --- |
| $D=0$ | [two_electron_unchanged_occupations.f90](two_electron_unchanged_occupations.f90)，一壳或两壳占据保持 |
| $D=2$ | [two_electron_single_transfer.f90](two_electron_single_transfer.f90)，净转移一个电子 |
| $D=4,\ G=40$ | [two_electron_pair_transfer.f90](two_electron_pair_transfer.f90)，两个电子成对转移 |
| $D=4,\ G=22$ | [two_electron_three_shell_transfer.f90](two_electron_three_shell_transfer.f90)，三壳层间分出或汇入一对电子 |
| $D=4,\ G=4$ | [two_electron_four_shell_transfer.f90](two_electron_four_shell_transfer.f90)，按四壳层位置选择配对 |

单壳占据差超过 2 或 $D>4$ 时无双电子贡献。开壳层路径可临时将核心轨道加入耦合链。CSF 对角元中启用闭壳层项时，其余核心项通过回调处理。例如 Coulomb 的 [`cord`](../core/closed_shell_mcp.f90) 对同轨道单极项使用 $N(N-1)/2$，不同轨道的直接单极项使用 $N_aN_b$。

### Coulomb 球张量因子

[coulomb_angular_factor.f90](coulomb_angular_factor.f90) 检查 $(l_a,l_c,k)$、$(l_b,l_d,k)$ 的三角关系与偶宇称。通过后返回

$$
\begin{aligned}
C_k(a,b,c,d)&=
\eta_4(l_c-2k-l_a+l_d-l_b)\,
\mathcal C(|\kappa_a|,k,|\kappa_c|)\,
\mathcal C(|\kappa_b|,k,|\kappa_d|),\\
\mathcal C(u,k,v)&=(-1)^u\sqrt{4uv}\ \mathrm{clrx}(u,k,v),
\qquad u,v>0 .
\end{aligned}
$$

$\mathcal C$ 对应 [`cre`](../core/reduced_angular_matrix_element.f90)，$\mathrm{clrx}$ 来自[相对论球谐张量的 $3j$ 计算](../core/relativistic_three_j.f90)。$C_k$ 只含角向量，径向部分在装配时与之配对。

### 子壳层乘积的变换与相位

设 $M_\lambda$ 是对应分支的 `gg*` 或 `ww1` 子壳层乘积，$R_\lambda$ 为该分支选用的重耦合因子。代码为不同轨道配对重复使用以下两种形式：

$$
\begin{aligned}
\mathcal D(k)&=\frac{M_kR_k}{\sqrt{[k]}},\\
\mathcal U(k;p,q;r,s;\theta)&=
\sum_{\lambda}\eta_4(\theta)
\sqrt{[\lambda]}
\begin{Bmatrix}j_p&j_q&k\\j_r&j_s&\lambda\end{Bmatrix}
M_\lambda R_\lambda .
\end{aligned}
$$

$\lambda$ 遍历共同的三角区间，再通过 $6j$ 三角与分支中的相位奇偶选择。准备好的 $M_\lambda$、$R_\lambda$ 按源码的 `rank_table_position`（$\lambda+1$）保存。为读出各分支的正负号，下文用 $t_s=2j_s$，

$$
\sigma_l=\eta_4\!\left(\sum_s l_s d_s\right),\qquad
s_{pq}=\eta_2(N_{\mathrm{span}}(p,q)+1).
$$

$N_{\mathrm{span}}(p,q)$ 是从两个活跃位置中较早者到较晚者之前的 bra 占据数之和。[subshell_permutation_sign.f90](subshell_permutation_sign.f90) 用四个壳层位置的逆序数给出 $s_{\mathrm{perm}}=(-1)^{\#\{p<q:i_p>i_q\}}$。

下表中的 $V^k$ 指传给 `speak` 的角向系数。写成 $V^k=\sigma_l C_k B_k$，其中 $C_k$ 按源码相应的四轨道顺序取值。

| 路径 | $M_\lambda$ 与 $R_\lambda$ | $B_k$ |
| --- | --- | --- |
| 两壳层占据相同，直接项 | `gg1122` 与 $R_2$ | $\mathcal D(k)$ |
| 两壳层占据相同，交叉项 | 同一组 `gg1122` 与 $R_2$ | $\mathcal U(k;1,2;2,1;0)$ |
| 两壳层成对转移 | `gg1122` 与 $R_2$ | $-\tfrac12\,\mathcal U(k;1,2;2,1;t_1+t_2+2k+2\lambda)$ |
| [第一壳层有三个张量](single_transfer_triple_on_first_shell.f90) | `gg1222` 按 $(2,1)$ 取壳，$R_2$ 的秩为 $j_2$ | $-s_{12}\mathcal U(k;2,1;1,1;2t_1+2k+2\lambda)$ |
| [第二壳层有三个张量](single_transfer_triple_on_second_shell.f90) | `gg1112` 按 $(2,1)$ 取壳，$R_2$ 的秩为 $j_1$ | $-s_{12}\mathcal U(k;2,2;1,2;2t_2+2k+2\lambda)$ |
| [另有一个占据不变的壳层](single_transfer_with_spectator_shell.f90) | `gg1233` 与按 $(2,1,3)$ 计算的置换 $R_3$ | 直接路径为 $s_{12}s_{\mathrm{perm}}\mathcal D(k)$，交叉路径为 $s_{12}s_{\mathrm{perm}}\mathcal U(k;1,3;3,2;0)$ |
| [三壳层分出或汇入一对](three_shell_pair_split_or_join.f90) | `gg1233` 与按 $(1,2,3)$ 计算的置换 $R_3$ | $-s_{12}s_{\mathrm{perm}}\mathcal U(k;3,1;2,3;\theta)$。两种排序的 $\theta$ 分别为 $t_3+t_1+2\lambda+2k$、$t_2+t_3+2\lambda+2k$ |

表中 $s_{\mathrm{perm}}$ 在单转移有旁观壳层时按 $(2,1,3,3)$ 计算，在三壳层分出或汇入一对时按 $(3,1,2,3)$ 计算，其中数字指例程中的壳层参数位置。双电子入口只向这些分支传入相容的占据。三张量分支还分别要求 $t_2-\lambda+1$ 或 $t_1-\lambda+1$ 为偶数；三壳层分出或汇入一对时根据排序在两者中选择。

[two_electron_unchanged_occupations.f90](two_electron_unchanged_occupations.f90) 的单壳层路径改用 $M_\lambda=(WW)^0_\lambda$。以该壳层角动量 $j$ 为例，

$$
V^k=
-\frac{C_k(a,a,a,a)}{2\sqrt{[J_{L,a}]}}
\sum_\lambda
\eta_2(2j+\lambda+k)\sqrt{[\lambda]}
\begin{Bmatrix}j&j&k\\j&j&\lambda\end{Bmatrix}
(WW)^0_\lambda .
$$

### 四壳层配对的计算式

四壳层三个求值文件都使用 $M_\lambda=A_1A_2A_3A_4$ 和 $R_4$，并取公共因子 $s_{12}s_{34}$。下表给出 $B_k/(s_{12}s_{34})$。设 $b$ 为 `orbital_order_case`，$(\theta_1,\theta_2)$ 分别对应 $b=1,2$。

| 文件 | 直接路径 | 交叉路径 |
| --- | --- | --- |
| [four_shell_transfer_12_34.f90](four_shell_transfer_12_34.f90) | $-\mathcal U(k;1,3;4,2;\theta)$，$\theta_1=t_2+t_3+2\lambda+2k$，$\theta_2=t_1+t_4-2\lambda+2k$ | $\mathcal U(k;1,4;3,2;\theta)$，$\theta_1=t_2+t_3+2t_4+2k$，$\theta_2=t_1+2t_2+t_4+4\lambda+2k$ |
| [four_shell_transfer_13_24.f90](four_shell_transfer_13_24.f90) | $\mathcal D(k)$，$b=2$ 时再乘 $\eta_4(t_1+t_2+t_3+t_4+4k)$ | $-\mathcal U(k;1,4;3,2;\theta)$，$\theta_1=2t_3-4k-4\lambda$，$\theta_2=t_1+t_2+t_3-t_4-4k$ |
| [four_shell_transfer_14_23.f90](four_shell_transfer_14_23.f90) | $\mathcal D(k)\eta_4(\theta)$，$\theta_1=t_4+t_3-2k+2$，$\theta_2=t_1+t_2-2k+2$ | $\mathcal U(k;1,3;4,2;\theta)$，$\theta_1=t_3-t_4-4k+2\lambda$，$\theta_2=t_1-t_2-4k-2\lambda$ |

以上四壳层直接与交叉路径的相位已分开列出。[分派入口](two_electron_four_shell_transfer.f90) 从轨道位置选择配对与 $b$。各求值例程按 $b$ 决定发给 `coulom` 和 `speak` 的轨道顺序。

### Breit 径向阶数和系数累加

[breit_radial_rank_ranges.f90](breit_radial_rank_ranges.f90) 对四个轨道采用直接配对 $(1,3),(2,4)$ 或交换配对 $(1,4),(2,3)$。对任一配对 $(a,b)$，令 $p_{ab}=1$ 当 $\kappa_a\kappa_b<0$，否则为 0，计算

$$
\begin{aligned}
L_{ab}&=\tfrac12|[j_a]-[j_b]|-1+p_{ab},\\
U_{ab}&=\tfrac12([j_a]+[j_b])+p_{ab}.
\end{aligned}
$$

若 $L_{ab}=-1$，代码取 $L_{ab}=1$。两个下界奇偶不同则拒绝该配对。令 $L$ 为两个下界的较大值，$U$ 为两个上界的较小值，则准许的 $\nu$ 为 $L,L+2,\ldots\le U$。代码返回起点 $L$ 和整数项数 $1+\operatorname{trunc}((U-L)/2)$，对应 `first_direct_nu`、`direct_nu_count` 或交换路径的两项。同侧两个轨道标签相同时跳过交换范围，其余轨道相等关系用于选择积分子类。

对固定多极阶数 $k$，普通多壳层路径对允许的 $\nu\in\{k-1,k,k+1\}$ 调用 [`calculate_breit_radial_integral_coefficients`](../breit_angular_support/breit_radial_integral_coefficients.f90)。设 $A_k^{D/E}$ 为对应的重耦合结果（例如四壳层表中的 $B_k$），$s_\mu^{D/E}(\nu,k)$ 为该过程返回的第 $\mu$ 个张量因子。`direct_breit_coefficients` 和 `exchange_breit_coefficients` 实现

$$
B_{\mu,\nu}^{D/E}=\sum_{k:\,|\nu-k|\le1}
A_k^{D/E}\,s_\mu^{D/E}(\nu,k).
$$

累加前按各分支检查轨道三角、径向阶数奇偶与下标范围。数组列号为 $1+(\nu-L)/2$。例如四壳层直接路径的 $\mu=1,\ldots,4$ 按径向阶数 $\nu$ 和积分类型 1 传给 `collect_breit_angular_coefficient`，$\mu=5,\ldots,12$ 按 $\nu+1$ 和类型 2 传给 `collect_breit_angular_coefficient`，各槽带有自己的轨道排列。单壳层路径直接发出 $-s_1 A_k/(2\sqrt{[J_{L,a}]})$；有三个张量落在同一壳层的单转移分支直接发出 $-s_1 A_k$。具体分支的排列和索引见[单壳层与同占据](two_electron_unchanged_occupations.f90)、[第一壳层三张量](single_transfer_triple_on_first_shell.f90)、[第二壳层三张量](single_transfer_triple_on_second_shell.f90)及[四壳层积分累加](four_shell_transfer_12_34.f90)。

## 来源

角向实现源自 GRASP 的 Gaigalas 例程。Racah 求值实现的原始注释记录了 N. S. Scott 与 C. Froese Fischer 的贡献。GRASP 许可证见[仓库 LICENSE](../../../LICENSE)。

G. Gaigalas 于 2009-11-22 记录了 GRASP2K 中两处需要比较 CSF 编号的修正。仅 seniority 不同的两个 CSF 可以有相同的占据数，仍应按非对角元处理。现行 `calculate_one_electron_scalar_coefficients` 在这种情况下将同壳层标量系数置零，`calculate_two_electron_angular_coefficients` 只在 CSF 对角元调用末段的闭壳层回调。 [单电子标量入口](one_electron_scalar_coefficients.f90) · [双电子角向入口](two_electron_angular_coefficients.f90)
