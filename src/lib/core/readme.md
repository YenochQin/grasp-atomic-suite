# Core 公共基础库

`core` 为 GRASP 的应用程序提供 CSF 基组输入、径向数值计算、角动量运算和通用 I/O，构建产物为 `libcore.a`。本文用公式说明计算例程的实际实现。其他例程的入口行为见各文件开头的注释。文件、模块和变量的旧名见[命名迁移](../../../docs/implemented/core_naming.md)与[变量迁移表](../../../docs/implemented/core_variable_migration.md)。

下文的 $P_a(r), Q_a(r)$ 是轨道 $a$ 的大小径向分量，$\kappa_a$ 是其相对论角量子数。径向量使用原子单位，$c_0$ 表示原子单位下的光速，避免与 Fermi 核半密度半径 $c_F$ 混淆。

## 径向网格、导数和求积

### 变换网格

[`radial_grid_setup.f90`](radial_grid_setup.f90) 的 `radgrd` 在等间隔变换坐标 $t_i=h(i-1)$ 上生成半径 $r_i$。$h$ 对应 `transformed_grid_step`，$h_p$ 对应 `asymptotic_grid_step`，$r_0$ 对应 `radial_grid_scale`。

当 $h_p=0$ 时，

$$
r(t)=r_0(e^t-1), \qquad \frac{dr}{dt}=r_0e^t.
$$

当 $h_p\ne0$ 时，按网格顺序求解隐式方程。令 $\eta=h/h_p$，则

$$
F(r;t)=\ln(1+r/r_0)+\eta r-t=0,\qquad
\frac{dr}{dt}=\frac{r+r_0}{1+\eta(r+r_0)}.
$$

每一步 Newton 修正为 $\Delta r=-F(r;t)(dr/dt)$，满足 $\lvert\Delta r/r_{\rm new}\rvert<10^3\,\mathtt{precis}$ 时接受。例程保存 $r_i$、$r'_i$ 和 $r'_i/r_i$，其中 $r'_i=dr/dt|_{t_i}$。原点的比值在数组中置零，活动区之后还生成 10 个网格点供工作数组使用。

### 差分权重与步长缩放的导数

[`radial_stencil_weights.f90`](radial_stencil_weights.f90) 和 [`radial_wavefunction_derivatives.f90`](radial_wavefunction_derivatives.f90) 使用 13 个等间隔节点。对单位间隔节点 $u_m=m-1$，定义 Lagrange 基函数与导数权重

$$
L_m(u)=\prod_{\substack{\nu=1\\\nu\ne m}}^{13}
   \frac{u-u_\nu}{u_m-u_\nu},
\qquad A_{jm}=L'_m(u_j).
$$

源码以公共分母 $479001600$ 保存 $A$。`dpbdt` 对 $F=P,Q$ 分别计算

$$
h\frac{dF}{dt}(t_i)\simeq \sum_{m=1}^{13} A_{jm}F(t_{i_m}).
$$

前六点用最前 13 个采样点及 $j=i$。最后六点用末尾 13 个采样点及 $j=i-N+13$。中间用第 7 行的反对称模板，等价于
$\sum_{s=1}^6 A_{7,7-s}[F(t_{i-s})-F(t_{i+s})]$。

同一个初始化例程还为 5 点和 6 点闭合 Newton–Cotes 求积生成权重。若 $L_m^{(n)}$ 是 $n$ 点单位间隔的 Lagrange 基函数，则从第 1 点积分到第 $j$ 点的权重为

$$
C^{(n)}_{mj}=h\int_0^{j-1}L_m^{(n)}(u)\,du.
$$

5 点部分区间使用 $j=2,3,4$，分母为 720。Hartree 内积分的原点被积函数为零，因此这个工作表只存储权重的第 2 到第 5 行。6 点使用 $j=2,\ldots,6$，分母为 1440。完整 5 点面板的权重为

$$
(c_1,c_2,c_3,c_2,c_1)=\frac{h}{45}(14,64,24,64,14).
$$

相邻面板共享端点，合并端点权重为 $c_4=2c_1$。径向步进另有一组起始权重 $(1771,9235h,5890h,4610h,35h,59h)/5760$。

### 原点拟合与复合求积

[`radial_quadrature.f90`](radial_quadrature.f90) 的 `quad` 接收 $T_i=r'_i f(r_i)$。从 $i=2$ 起寻找首对同号非零样本，以 $f(r)\simeq Ar^\sigma$ 拟合小半径区间。

$$
\sigma=\frac{\ln\!\big[(T_{i+1}/r'_{i+1})/(T_i/r'_i)\big]}
              {\ln(r_{i+1}/r_i)}, \qquad
I_{0,i}\simeq\frac{r_i (T_i/r'_i)}{\sigma+1}.
$$

余下区间对 $T_i$ 用上述 5 点面板累加。例程把最后有效点后的三个工作数组元素置零，以便覆盖最后一个面板。若最后有效点恰与面板端点对齐，减去一次 $c_1 T_{\rm last}$ 以校正端点权重。无可用拟合点对时返回零。

## 轨道插值、正交化与节点

[`radial_wavefunction_records.f90`](radial_wavefunction_records.f90) 按 $(n,\kappa)$ 匹配轨道记录。载入前设置原点指数，点核为 $\gamma=\sqrt{\kappa^2-(Z/c_0)^2}$，有限核为 $\gamma=|\kappa|$。从文件读出的 $P,Q$ 先插值到当前网格，再对整套已载入轨道执行正交化。

### Aitken 插值

[`aitken_interpolation.f90`](aitken_interpolation.f90) 和 [`radial_orbital_interpolation.f90`](radial_orbital_interpolation.f90) 在靠近目标 $x$ 的搜索窗口中逐次选取尚未使用的节点。设选中节点为 $(x_i,y_i)$，令 $\delta_i=x_i-x$。三角表从 $T_{i,1}=y_i$ 开始，并用

$$
T_{i,k}=
\frac{T_{k-1,k-1}\delta_i-T_{i,k-1}\delta_{k-1}}{x_i-x_{k-1}},
\quad 2\le k\le i,\qquad
y^{(i)}(x)=T_{i,i}.
$$

`interp` 最多用 11 点，从相邻估计中选相对变化最小的估计。`intrpq` 对 $P,Q$ 共用节点，最多用 13 点。两个新分量非零时，只有各自的 $\lvert(y^{(i)}-y^{(i-1)})/y^{(i)}\rvert$ 都低于容差才提前接受。某个估计为零时走节点上限分支。它将超出有效范围的尾部清零，并用下式归一化，返回归一化前的范数平方。

### 重叠积分与 Gram–Schmidt

[`weighted_radial_overlap.f90`](weighted_radial_overlap.f90) 的 `rint` 对两轨道共同的有效径向范围求

$$
S_{ab}^{(k)}=\int_0^\infty r^k
   [P_a(r)P_b(r)+Q_a(r)Q_b(r)]\,dr .
$$

[`radial_orbital_orthogonalization.f90`](radial_orbital_orthogonalization.f90) 依次处理轨道。令 $u_a=(P_a,Q_a)$，对每个已处理且 $\kappa_b=\kappa_a$ 的轨道，用当前剩余分量计算 $S_{ab}^{(0)}$，再更新

$$
u_a\leftarrow u_a-S_{ab}^{(0)}u_b,\qquad
u_a\leftarrow\frac{u_a}{\sqrt{S_{aa}^{(0)}}}
\quad\text{（做过投影时，在所有投影后归一化）}.
$$

原点系数跟随同一投影和归一化更新。归一化后，从有效区末端向内清除 $\lvert P_a\rvert<0.01\,\mathtt{accy}$ 的尾点。

### 径向节点计数

[`radial_node_count.f90`](radial_node_count.f90) 在跨号处分段，并记录每段绝对振幅最大的样本。以 Fortran `SIGN` 实现符号函数 $s$，分段和后续节点计数的判据都是 $\lvert s(F_i)+s(F_{i-1})\rvert\le\mathtt{accy}$。若分段极值的幅度小于

$$
A_{\min}=\mathtt{thresh}\,\max_{2\le i\le i_{\rm last}}\lvert F_i\rvert ,
$$

则移除该极值。在保留的相邻极值之间重新计数，并另返回首个保留极值的符号。

## Dirac–Coulomb 初猜和径向方程

### 点 Coulomb 场的解析束缚轨道

[`dirac_coulomb_radial_orbital.f90`](dirac_coulomb_radial_orbital.f90) 的 `dcbsrw` 使用有效电荷 $Z_{\rm eff}$、主量子数 $n$ 和 $\kappa$。令

$$
n_r=n-|\kappa|,\quad z=Z_{\rm eff}/c_0,\quad
\gamma=\sqrt{\kappa^2-z^2},\quad
\mathcal N=\sqrt{n^2-2n_r(|\kappa|-\gamma)} .
$$

返回的束缚能 $E_b$ 以连续谱阈值为零点并取正值。

$$
\varepsilon=\left[1+\left(\frac{z}{n_r+\gamma}\right)^2\right]^{-1/2},
\qquad E_b=c_0^2(1-\varepsilon),
\qquad \rho=\frac{2 Z_{\rm eff}r}{\mathcal N}.
$$

设 $(a)_m$ 为升阶乘，$M(-q,b,\rho)=\sum_{m=0}^q(-q)_m\rho^m/[(b)_m\,m!]$。两个终止级数及公共径向因子为

$$
F_1=n_r M(1-n_r,2\gamma+1,\rho),\quad
F_2=(\mathcal N-\kappa) M(-n_r,2\gamma+1,\rho),\quad
W=e^{-\rho/2}\rho^\gamma .
$$

$n_r=0$ 时 $F_1=0$。按源码的相位约定，

$$
\begin{aligned}
A&=\sigma_\kappa
  \frac{\sqrt{\Gamma(2\gamma+1+n_r)}}{\Gamma(2\gamma+1)\sqrt{n_r!}}
  \sqrt{\frac{Z_{\rm eff}}{2\mathcal N^2(\mathcal N-\kappa)}},
&\sigma_\kappa&=\begin{cases}+1&\kappa>0,\\-1&\kappa<0.\end{cases}\\
P&=A\sqrt{1+\varepsilon}\,W(F_1-F_2),
&Q&=A\sqrt{1-\varepsilon}\,W(F_1+F_2).
\end{aligned}
$$

原点大分量系数为 $A\sqrt{1+\varepsilon}(2Z_{\rm eff}/\mathcal N)^\gamma[n_r-(\mathcal N-\kappa)]$。末尾满足 $|P_i|<0.1\,\mathtt{accy}$ 或 $|r_i|>10^{50}$ 的采样点被清零。

### 复 Gamma

[`complex_gamma.f90`](complex_gamma.f90) 将参数平移到级数适用区，再通过递推乘积恢复原参数。令 $w=z+m$，其中 $m$ 是使 $|w|\ge10$ 的平移步数，则

$$
\begin{aligned}
L(w) &\simeq (w-\tfrac12)\log w-w+\tfrac12\log(2\pi)
 +\sum_{j=1}^{7}\frac{B_{2j}}{2j(2j-1)\,w^{2j-1}},\\
\Gamma(z)&\simeq \exp\!\left(L(w)-\log\prod_{\ell=0}^{m-1}(z+\ell)\right).
\end{aligned}
$$

负实参数通过 $\Gamma(x)\Gamma(-x)=-\pi/[x\sin(\pi x)]$ 恢复，接近整数极点时停止。[`nonnegative_atan2.f90`](nonnegative_atan2.f90) 提供复对数所用的相角分支。例程最后把 $L$ 的实部指数化，再用虚部的正弦和余弦组成实、虚结果。

### 网格上的 Dirac 系数和匹配点

[`radial_dirac_coefficients.f90`](radial_dirac_coefficients.f90) 将 `yp` 表中的轨道势量记为 $Y_i$，将 `orbital_energy_parameters` 记为 $E$。对 $i\ge2$ 写入径向步进工作区

$$
T_{F,i}=-\frac{h r'_i}{2c_0}\left(2c_0^2-E+\frac{Y_i}{r_i}\right),
\qquad
T_{G,i}=-\frac{h r'_i}{2c_0}\left(E-\frac{Y_i}{r_i}\right).
$$

原点的两项置零。以 Fortran `SIGN` 找到首个 $T_{G,i}$ 与 $T_{G,1}$ 反号的点，具体判据是两符号之和的绝对值小于 `accy`。找不到时使用编译期最大径向点数的一半。

### 原点六点迭代

[`radial_origin_solution.f90`](radial_origin_solution.f90) 对 $P_i=r_i^\gamma p_i,\ Q_i=r_i^\gamma q_i$ 的缩放分量迭代。定义 $g_\pm=1-\gamma\pm\kappa$，并令

$$
D_{P,i}=-\frac{r'_i}{c_0}\big[(2c_0^2-E)r_i+Y_i\big],
\qquad
D_{Q,i}=-\frac{r'_i}{c_0}\big[E r_i-Y_i\big].
$$

第 $j=2,\ldots,6$ 点使用 6 点部分区间权重更新

$$
\begin{aligned}
p_j^{\rm new}
 &=\frac1{r_j}\sum_{i=1}^{6} C^{(6)}_{ij}
        (g_- r'_i p_i-D_{P,i}q_i+S_{P,i}),\\
q_j^{\rm new}
 &=\frac1{r_j}\sum_{i=1}^{6} C^{(6)}_{ij}
        (g_+ r'_i q_i-D_{Q,i}p_i+S_{Q,i}).
\end{aligned}
$$

令 $s$ 为传入的原点斜率参数，$z=Z/c_0$。点核齐次及非齐次模式均用 $p_1=s$。$\kappa<0$ 时 $q_1=-sz/(\gamma-\kappa)$，否则 $q_1=s(\gamma+\kappa)/z$。有限核在 $\kappa<0$ 时改用 $(p_1,q_1)=(s,0)$，否则为 $(0,s(\gamma+\kappa)/z)$。变分模式的两个原点系数均为零。

齐次模式 $S=0$。非齐次模式在 $i\ge2$ 用
$(S_{P,i},S_{Q,i})=-r'_i r_i^{-\gamma}(X_{P,i},X_{Q,i})$，其中 $X_P,X_Q$ 是交换源。原点的源由直接项系数 $d_\ell$ 加权，其中 $(p_{1,\ell},q_{1,\ell})$ 用相应轨道的原点系数和当前方程的 $\gamma,\kappa,z$ 按上面的核模型规则构造。令 $\chi=1$ 对应点核，$\chi=0$ 对应有限核，则

$$
(S_{P,1},S_{Q,1})
 = r'_1 \sum_\ell d_\ell\,
   (g_-p_{1,\ell}+\chi zq_{1,\ell},\ g_+q_{1,\ell}-\chi zp_{1,\ell}).
$$

变分模式在 $i\ge2$ 取 $r'_i r_i^{1-\gamma}(-Q_i^{\rm old},P_i^{\rm old})/c_0$，原点源为零。用 $r_6^\gamma\max_j(|\Delta p_j|,|\Delta q_j|)$ 与 `accy` 比较，最多迭代 36 次。超过上限时以初始原点系数填充第 2 到第 6 点，随后乘回 $r^\gamma$。输出的第 1 个物理采样点置零，小分量原点系数单独返回。

## 点核、Fermi 核及核极化

### 电子所见核势

[`nuclear_coulomb_potential.f90`](nuclear_coulomb_potential.f90) 保存 $\mathcal Z(r)=-rV(r)$。点核直接令 $\mathcal Z=Z$。对于半密度半径 $c_F$、弥散度 $a$ 的 Fermi 核，令

$$
u=\frac{a}{c_F},\quad
x=\frac{r}{c_F},\quad
y=\frac{r-c_F}{a},\quad
S_k(t)=\sum_{m=1}^{\infty}\frac{(-1)^m e^{mt}}{m^k},\quad
d=-6u^3 S_3(-1/u),\quad \mathcal A=1+\pi^2u^2+d .
$$

[`alternating_exponential_series.f90`](alternating_exponential_series.f90) 同步累加 $S_2,S_3$。当新旧 $S_2$ 的浮点绝对值相等时停止。实际求势分内外两段。

$$
\begin{aligned}
\mathcal Z(r\le c_F)
 &=\frac{Z}{\mathcal A}
   \left[d+6u^3S_3(y)+x\!\left(\frac32+\frac{\pi^2u^2}{2}
                          -3u^2S_2(y)-\frac{x^2}{2}\right)\right],\\
\mathcal Z(r>c_F)
 &= Z\left[1+\frac{3u^2}{\mathcal A}
             \left(x S_2(-y)+2u S_3(-y)\right)\right].
\end{aligned}
$$

### 核极化修正

[`nuclear_polarization_potential.f90`](nuclear_polarization_potential.f90) 按开启的通道求和。

$$
V_{\rm E1}(r)=-\frac{K_1}{r^4+b_1^4},\qquad
V_{\rm E2}(r)=-\frac{K_2}{r^6+b_2^6},\qquad
\mathcal Z_{\rm total}(r)=\mathcal Z_{\rm Coul}(r)-r\,V_{\rm NP}(r).
$$

其中 $V_{\rm NP}$ 是开启通道之和。为避免直接计算大幂，令 $s=\max(r,b)$，把 $r^p+b^p$ 写为 $s^p[(r/s)^p+(b/s)^p]$。`channel` 用二进制尾数和指数分别表示 $K/s^p$，并检查结果是否可表示。原点合成项为零。

## 单电子与双电子径向积分

### 单电子 Dirac 径向积分

[`one_electron_radial_integral.f90`](one_electron_radial_integral.f90) 要求 $\kappa_a=\kappa_b$。记 $D_{ab}=P_aP_b+Q_aQ_b$，则三项动能与势能项在连续形式下为

$$
\begin{aligned}
I_{ab}={}&
 c_0\int_0^\infty\left(Q_a\frac{dP_b}{dr}-P_a\frac{dQ_b}{dr}\right)\,dr
 -2c_0^2\int_0^\infty Q_aQ_b\,dr\\
 &+c_0\kappa_b\int_0^\infty\frac{P_aQ_b+Q_aP_b}{r}\,dr
 -\mathbf 1_{\mathrm{potential\_mode}=0}
       \int_0^\infty \frac{\mathcal Z(r)}{r}D_{ab}(r)\,dr .
\end{aligned}
$$

代码第一项用 `dpbdt` 给出的 $h\,dP/dt,h\,dQ/dt$ 写入求积数组并除以 $h$，后三项按 $r'_i$ 或 $r'_i/r_i$ 加权。势能项读取已装入的 $\mathcal Z$，可包含启用的核极化修正。

### Hartree 积分和 Slater 积分

[`hartree_yz_functions.f90`](hartree_yz_functions.f90) 以 $D_{bd}(s)=P_b(s)P_d(s)+Q_b(s)Q_d(s)$ 构造

$$
\begin{aligned}
Z^k_{bd}(r) &= r^{-k}\int_0^r s^k D_{bd}(s)\,ds,\\
W^k_{bd}(r) &= r^{k+1}\int_r^\infty s^{-k-1}D_{bd}(s)\,ds,\\
Y^k_{bd}(r) &= Z^k_{bd}(r)+W^k_{bd}(r).
\end{aligned}
$$

内积分从原点向外、外积分从尾部向内，均以 5 点面板递推。$k=0$ 时以 $\delta_{bd}$ 设定内积分的尾部极限。$k>0$ 时用最后的内积分乘以 $r^{-k}$ 延拓。网格上的径向幂按张量阶数缓存。

[`slater_radial_integral.f90`](slater_radial_integral.f90) 将 $Y^k_{bd}$ 与另一对轨道积相乘并积分。

$$
R^k(ab;cd)=\int_0^\infty\frac{D_{ac}(r)\,Y^k_{bd}(r)}{r}\,dr
=\int_0^\infty\!\!\int_0^\infty
 D_{ac}(r)D_{bd}(s)\frac{\min(r,s)^k}{\max(r,s)^{k+1}}\,ds\,dr .
$$

## 角动量、CSF 与 Coulomb 角系数

### 相对论 3-j 与约化矩阵元

[`log_factorial_table.f90`](log_factorial_table.f90) 提供 $g(s)=\ln[(s-1)!]$。[`relativistic_three_j.f90`](relativistic_three_j.f90) 计算

$$
w_k(a,b)=
\begin{pmatrix} j_a&k&j_b\\ \frac12&0&-\frac12 \end{pmatrix},
\qquad j_i=|\kappa_i|-\tfrac12 .
$$

令 $a=|\kappa_a|,\ b=|\kappa_b|,\ A=a+b,\ d=a-b,\ S=A+k$。若 $k<|d|$ 或 $k>A-1$，返回零。在允许区间，公共对数项为

$$
L_0=\tfrac12[-\ln(ab)+g(A-k)+g(d+k+1)+g(-d+k+1)-g(S+1)].
$$

当 $S$ 为偶数时取 $(\delta,\eta,p)=(2,0,3S/2)$，否则取 $(1,1,(3S-1)/2)$。源码的结果等价于

$$
w_k=(-1)^p \exp\!\left[
 L_0+g\!\left(\frac{S+\delta}{2}\right)-g\!\left(\frac{A-k+\eta}{2}\right)
    -g\!\left(\frac{d+k+\delta}{2}\right)-g\!\left(\frac{-d+k+\delta}{2}\right)
 \right].
$$

[`reduced_angular_matrix_element.f90`](reduced_angular_matrix_element.f90) 乘上简并度和相位，返回
$\mathrm{CRE}=(-1)^a\sqrt{(2j_a+1)(2j_b+1)}\,w_k(a,b)$。

### Coulomb 允许阶数与闭壳层配对

[`coulomb_tensor_rank_range.f90`](coulomb_tensor_rank_range.f90) 对一对轨道 $i,j$ 以 $s_i=2j_i+1$ 定义二倍阶数上下界

$$
2L_{ij}=|s_i-s_j|+2\mathbf 1_{\kappa_i\kappa_j<0},
\qquad
2U_{ij}=s_i+s_j-2-2\mathbf 1_{\kappa_i\kappa_j>0}.
$$

直接项比较 $(1,3),(2,4)$，交换项比较 $(1,4),(2,3)$。当两组下界的奇偶性相同时，首阶取 $k_{\min}=\max(L_1,L_2)$，末阶取 $k_{\max}=\min(U_1,U_2)$，步长为 2。程序返回的计数为 Fortran 整数式 $(k_{\max}-k_{\min})/2+1$。任一侧的两个轨道相同时跳过交换分支。

[`closed_shell_mcp.f90`](closed_shell_mcp.f90) 在同一轨道 $a=b$ 的 $k=0$ 直接项使用 $q_a(q_a-1)/2$，高阶偶数项使用 $-\lfloor q_a^2/2\rfloor\,w_k(a,a)^2$。两轨道不同的 $k=0$ 直接项使用 $q_aq_b$，交换项使用 $-q_aq_b\,w_k(a,b)^2$，其阶数从 $\big||\kappa_a|-|\kappa_b|\big|+\mathbf 1_{\kappa_a\kappa_b<0}$ 起每次加 2，止于 $|\kappa_a|+|\kappa_b|-1$。其中 $q_a,q_b$ 为占据数。交换项绝对值低于 $10^{-10}$ 时不写入缓冲。

### CSF 解码、筛选和重耦合表

[`csf_descriptor_reader.f90`](csf_descriptor_reader.f90) 接收每个活动（peel）子壳层的四元组 $(q,2J,v,2K)$，行末保存总 $2J$ 和宇称。完整行宽为 $4n_{\rm peel}+2$。对单电子 $\kappa$ 和各已占据子壳层，检查或推导

$$
l=\begin{cases}\kappa&\kappa>0,\\-\kappa-1&\kappa<0,\end{cases}
\qquad 2j=2|\kappa|-1,\qquad
N_e=\sum_i q_i,\qquad
\pi=(-1)^{\sum_i q_i l_i}.
$$

闭核的 $q_i=2j_i+1$，其偶宇称因子为 1。沿已占据子壳层连接中间角动量，每步须满足
$\lvert2K_{\rm prev}-2J_i\rvert\le2K_{\rm new}\le2K_{\rm prev}+2J_i$，且 $2K_{\rm prev}+2J_i-2K_{\rm new}$ 为偶数。总对称性以 $\pi(2J+1)$ 写入有符号字节，由 [`csf_quantum_queries.f90`](csf_quantum_queries.f90) 分别取符号和绝对值。

[`one_electron_csf_filter.f90`](one_electron_csf_filter.f90) 对 $\Delta q_i=q_i^{\rm bra}-q_i^{\rm ket}$ 要求非零差数目为 0 或 2，且 $|\Delta q_i|\le1$。[`two_electron_csf_filter.f90`](two_electron_csf_filter.f90) 要求非零差数目不超过 4，且 $|\Delta q_i|\le2$。[`angular_momentum_triangle.f90`](angular_momentum_triangle.f90) 对传入的 $s_i=2J_i+1$ 检查 $|s_2-s_3|+1\le s_1\le s_2+s_3-1$。

[`subshell_parent_table.f90`](subshell_parent_table.f90) 用 $m=\min(q-1,s-(q-1))+1$ 选择电子或空穴的表项，$s=2j+1$。$m=1$ 时用第 1 行，否则行号为 $s(s-2)/8+m$。同一侧两次作用于同一子壳层时，先将第一个占据取为第二个占据减 1。[`csf_pair_quantum_numbers.f90`](csf_pair_quantum_numbers.f90) 与 [`angular_recoupling_tables.f90`](angular_recoupling_tables.f90) 用这些占据和耦合码组织左右态的耦合树，供外部重耦合计算使用。

## 稀疏矩阵运算与初始本征基

[`sparse_symmetric_matvec.f90`](sparse_symmetric_matvec.f90) 对每个输入向量计算 $y=Hx$。存储三角的一列 $j$ 末尾是对角元。对该列的每个非对角存储项 $H_{ij}$，[`sparse_dot_and_scatter.f90`](sparse_dot_and_scatter.f90) 同时计算

$$
y_j\mathrel{+}=H_{ij}x_i,\qquad
y_i\mathrel{+}=H_{ij}x_j .
$$

[`sparse_eigensolver_initial_basis.f90`](sparse_eigensolver_initial_basis.f90) 选取阶数 $n_s=\min(n_{\max},N_{\rm CSF})$ 的左上子矩阵。上三角的行 $i$、列 $j$ 打包到一基索引 $i+j(j-1)/2$，交由 `DSPEVX` 求最低的指定数量本征对。求得的子空间向量在其余 CSF 分量处补零，能量附在向量区之后。

## 物理常数换算

[`physical_constants.f90`](physical_constants.f90) 用两倍 Rydberg 值给出 Hartree 到 eV 及 $\mathrm{cm}^{-1}$ 的换算。由 Bohr 半径的厘米值 $a_0$ 得 $1\,\mathrm{fm}=10^{-13}/a_0$ Bohr。Einstein 系数的换算因子在实现中为

$$
C_A=\frac{m_e}{\hbar}\left(\frac{e_{\rm esu}^{2}}{\hbar}\right)^2,\qquad
C_B=\frac{10a_0^3}{\hbar}C_A ,
$$

其中输入为文件中的 CGS 常数。

## 文本径向图的坐标

[`radial_ascii_plot.f90`](radial_ascii_plot.f90) 先对两分量应用输入尺度，
求参与绘图的合并最小值 $m$ 和最大值 $M$。有非零幅值跨度时，行缩放因子为
$A=54/|M-m|$。如果 $m<0<M$，零基线行 $i_0$ 取 $4+A|m|$ 的 Fortran 整数截断值，
否则 $i_0=1$。第 $i$ 个网格点上已缩放的分量 $F_i$ 对应

$$
\mathrm{row}(F_i)=i_0+\operatorname{NINT}(AF_i),
\qquad
\mathrm{column}(i)=\max(1,\operatorname{NINT}(131\,i/N)).
$$

$N$ 是活动网格点数。因此横轴按采样点编号排布。两个分量全零时直接报告无曲线。
