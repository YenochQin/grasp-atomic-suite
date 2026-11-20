# Davidson 本征求解器

`solve_davidson_eigenpairs` 使用分块 Davidson 迭代，求实对称矩阵 $A$ 谱端的若干本征对。调用方给出 $A$ 的对角元、初始向量，以及对一批列向量执行 $AV$ 的回调 `dvd_operator`。完整矩阵的存储由调用方负责。下列公式按 `2d4c721` 基础上的工作树源码核对（2026-11-20）。[入口源码](davidson_eigensolver.f90)

| 文件 | 过程 | 职责 |
| --- | --- | --- |
| [davidson_eigensolver.f90](davidson_eigensolver.f90) | `solve_davidson_eigenpairs` | 执行矩阵乘法回调和对角预条件 |
| [davidson_solver_control.f90](davidson_solver_control.f90) | `dvdson` | 选择谱端、划分工作区、调度并打包结果 |
| [davidson_sparse_initial_estimates.f90](davidson_sparse_initial_estimates.f90) | `iniest` | 提取稀疏对称性块的主子矩阵并生成初猜 |
| [davidson_initial_subspace.f90](davidson_initial_subspace.f90) | `initdvd` | 扩充初始子空间并构造首个投影 |
| [davidson_orthonormalization.f90](davidson_orthonormalization.f90) | `mgs_nrm` | 对新增列执行两遍改进 Gram–Schmidt |
| [davidson_projected_matrix.f90](davidson_projected_matrix.f90) | `adds` | 追加上三角打包投影矩阵的列 |
| [davidson_iteration.f90](davidson_iteration.f90) | `dvdrvr` | 求解投影本征问题，扩展或重启子空间并回写结果 |
| [davidson_convergence.f90](davidson_convergence.f90) | `tstsel` | 检查 Ritz 值变化和最新基列上的系数，排序待处理态 |
| [davidson_residual_expansion.f90](davidson_residual_expansion.f90) | `newvec` | 计算 Ritz 残差及预条件右端项 |
| [davidson_subspace_restart.f90](davidson_subspace_restart.f90) | `ovflow` | 保留谱端 Ritz 向量并压缩子空间 |
| [davidson_basis_transform.f90](davidson_basis_transform.f90) | `multbc` | 按行暂存结果并原地右乘系数矩阵 |

## 目标本征对与工作算子

设 $N$ 为 `matrix_order`，$L$ 为 `subspace_limit`，$k$ 为当前子空间列数。原矩阵的本征值按升序编号。设 $a$ 和 $b$ 分别为请求序号的最小值和最大值。若范围起点有效，`dvdson` 用闭区间 $[a,b]$ 生成目标序号。否则读取 `selected_pair_indices`，遇到首个非正项时停止。

为使投影求解始终从低端取值，设 $\sigma$ 为谱端符号，$j$ 为原本征对序号，$j'$ 为内部序号，$m$ 为 `extreme_pair_count`。数值 $m$ 计算从所选谱端到最远目标的跨度。

$$
(\sigma,j',m)=
\begin{cases}
(+1,j,b), & a+b-1\le N,\\
(-1,N-j+1,N-a+1), & a+b-1>N.
\end{cases}
\qquad
\widehat A=\sigma A .
$$

以下的 Ritz 值 $\theta_i$ 属于 $\widehat A$，索引 $i$ 使用内部顺序。回调始终计算原矩阵 $A$ 的作用，初始化和迭代在把结果写入工作投影前应用 $\sigma$。[目标选择与结果打包](davidson_solver_control.f90)

## 稀疏初猜与初始子空间

设指定对称性块的阶数为 $n_b$。`iniest` 使用该块的前 $h=\min(1000,n_b)$ 个基函数。设它们在全空间的序号为 $p_1,\dots,p_h$，嵌入矩阵 $E=[e_{p_1},\dots,e_{p_h}]$ 由全空间坐标单位向量组成。提取出的主子矩阵和回填的初猜满足

$$
A_s=E^{\mathsf T}AE,\qquad
A_s y_i=\mu_i y_i,\qquad
b_i=Ey_i .
$$

对稀疏输入的第 $p_j$ 列，`column_end_offsets` 给出列首尾界限。若列内相邻存储项的全局行号为 $r_t,r_{t+1}$，令 $d_t$ 为代码中的 `row_index_shift`，$\beta$ 为目标块号，则实际索引更新为

$$
d_1=p_j-j,\qquad
d_{t+1}=d_t+\sum_{r=r_t}^{r_{t+1}}\mathbf{1}_{\operatorname{block}(r)\ne\beta},
\qquad
i_t=r_t-d_t .
$$

该存储项写入一基下三角打包位置 $i_t+(j-1)(2h-j)/2$。`DSPEVX` 求小矩阵的最低 `initial_vector_count` 个本征对。输出工作区把回填后的全空间列向量依次堆叠，再接上对应的 $\mu_i$。[稀疏初猜源码](davidson_sparse_initial_estimates.f90)

调用方也可提供已有初猜。`initdvd` 计算这些初始列的 $Ab_j$，取上一批矩阵乘积中所需的前几列正交化，作为下一批初始基。`dvdson` 传入的目标初始列数为 $m+1$。令 $B_{\rm batch}$ 为上一批基的相应列，$\operatorname{MGS}_2(V;Q)$ 表示将 $V$ 相对已有列 $Q$ 做两遍正交归一化，则新增基列为

$$
B_{\rm new}=\operatorname{MGS}_2(A B_{\rm batch};Q).
$$

直接初始化时，零初猜对应 $b_1=\mathbf 1/\sqrt N$。公开入口 `solve_davidson_eigenpairs` 要求调用方提供至少一个初猜。初始化结束时，依次应用谱端符号并构造第一份投影。[初始子空间源码](davidson_initial_subspace.f90)

### 两遍改进 Gram–Schmidt

对于每个已有的单位参考列 $q_i$，逐次从新增列 $v_j$ 去掉投影。随后依次归一化每个新列，并用归一化后的新列更新其后各列。下面的更新在两遍中各执行一次，归一化要求被处理的列具有非零范数。

$$
v_j\leftarrow v_j-q_i(q_i^{\mathsf T}v_j),
\qquad
q_j=\frac{v_j}{\sqrt{v_j^{\mathsf T}v_j}} .
$$

已有列包括给定约束和当前子空间基。公开入口用零个约束，迭代使用固定两遍正交化。[正交化源码](davidson_orthonormalization.f90)

## 投影和收敛判定

令 $B=[b_1,\dots,b_k]$ 为当前正交基，$W=\widehat A B$ 为已缓存的算子乘积。追加基列时只补算对应的新投影列。按列排列的上三角打包下标从 1 开始。

$$
S=B^{\mathsf T}W,\qquad
S_{ij}=b_i^{\mathsf T}W_{:j},\qquad
\pi_U(i,j)=i+\frac{j(j-1)}2,\quad i\le j .
$$

`dvdrvr` 将打包投影复制给 `DSPEVX`，求指定数量的低端投影本征对。设 $c_i$ 为投影空间的系数，则全空间近似向量为

$$
S c_i=\theta_i c_i,\qquad x_i=Bc_i .
$$

[投影源码](davidson_projected_matrix.f90) · [主迭代源码](davidson_iteration.f90)

`tstsel` 先检查所有选中态的 Ritz 值变化 $\Delta_i$。若每个 $\Delta_i<\varepsilon_E$，则结束迭代。否则对每个目标态计算最近一批基列上的绝对系数最大值 $g_i$。此处 $s$ 是提交给 `tstsel` 的新增列数，取值至少为 1 时使用对应的最后 $s$ 列。代码在 $s=0$ 时也保留最后一个系数参与检查。

$$
\Delta_i=|\theta_i-\theta_i^{\mathrm{old}}|,
\qquad
g_i=\max_{\ell=k-\max(1,s)+1,\ldots,k}|(c_i)_\ell| .
$$

低于系数阈值 $\varepsilon_C$ 的态标记为收敛。其余态按 $g_i$ 降序送去计算残差。若没有待处理态，系数判据结束迭代。[收敛检查源码](davidson_convergence.f90)

## 残差修正与对角预条件

`newvec` 利用同一组 Ritz 系数计算原始残差及其二范数。

$$
r_i=Wc_i-\theta_i Bc_i,\qquad
\rho_i=\sqrt{r_i^{\mathsf T}r_i}.
$$

若 $\rho_i<\varepsilon_R$，该候选态跳过扩展。待处理态都达到这一残差阈值时，残差判据结束迭代。否则用相邻 Ritz 值间隔 $\gamma_i$ 与 $\rho_i$ 比较，决定沿 $x_i$ 方向修正残差的幅度 $\eta_i$。此幅度使用前一轮保存的 Ritz 值变化 $\Delta_i$。

$$
\gamma_i=
\begin{cases}
|\theta_2-\theta_1|, & i=1,\\
\min(|\theta_i-\theta_{i-1}|,|\theta_{i+1}-\theta_i|), & i>1,
\end{cases}
\qquad
\eta_i=
\begin{cases}
\sqrt{\Delta_i\,\rho_i}, & \gamma_i>\rho_i,\\
\min(\Delta_i,\sqrt{\Delta_i}), & \gamma_i\le\rho_i .
\end{cases}
$$

$$
u_i=r_i-\eta_i x_i .
$$

$u_i$ 是预条件右端项。`newvec` 对每轮新增列数同时施加 `expansion_block_size` 和剩余子空间容量 $L-k$ 两个上限。[残差扩展源码](davidson_residual_expansion.f90)

`dvdrvr` 传给对角预条件器的原谱坐标值是 $\lambda_i^{\rm prec}=\sigma(\theta_i-\delta_i)$，其中 $\delta_i$ 是 `preconditioner_shifts(i)`，在当前 `newvec` 中为每个新增候选赋值为零。令 $D$ 为以 $A_{jj}$ 为对角元的对角矩阵。`solve_davidson_eigenpairs` 对第 $j$ 个分量计算

$$
t_{ji}=D_{jj}-\lambda_i^{\rm prec},
\qquad
(z_i)_j=
\begin{cases}
(u_i)_j/t_{ji}, & |t_{ji}|>10^{-6},\\
10^6(u_i)_j, & |t_{ji}|\le 10^{-6}.
\end{cases}
$$

`solve_davidson_eigenpairs` 在每个向量的预条件步骤结束时把位移加回对角数组。返回迭代后，$z_i$ 相对已有基和约束正交归一化，再由回调计算新增列的算子乘积。[预条件源码](davidson_eigensolver.f90) · [迭代调度](davidson_iteration.f90)

## 子空间重启与结果

`dvdrvr` 每轮从投影问题请求 $\min(k,d)$ 个低端本征对，其中重启保留列数 $d$ 由子空间容量和所需谱端跨度决定。

$$
d=\min\!\left(\max\!\left(\left\lfloor\frac L2\right\rfloor,2m\right),L-3\right).
$$

当 $k=L$ 且继续迭代时，`ovflow` 取最近一次投影求解中前 $d$ 个 Ritz 系数 $C_d=[c_1,\ldots,c_d]$。程序同步压缩基和算子乘积，在新基中直接将投影置为对角矩阵，并把投影本征向量设为单位阵。

$$
B'=BC_d,\qquad
W'=WC_d,\qquad
S'=\operatorname{diag}(\theta_1,\dots,\theta_d),\qquad
C'=I_d.
$$

`multbc` 每次暂存一整行乘积后才覆盖原数组的前 $d$ 列。用于返回结果时，乘积列数改为 $m$。满足收敛判据或子空间达到全空间阶数时，主迭代进入结果回写。超过迭代上限时回写上轮已求得的近似，并返回状态码。原谱中的本征值及用变换后的向量和乘积重新计算的残差范数满足

$$
\lambda_i=\sigma\theta_i,
\qquad
\|A x_i-\lambda_i x_i\|_2
=\|\widehat A x_i-\theta_i x_i\|_2 .
$$

[重启源码](davidson_subspace_restart.f90) · [基变换源码](davidson_basis_transform.f90)

## 反向通信与结果布局

| `operation_code` | 请求 |
| --- | --- |
| `0` | 开始或结束求解 |
| `1` | 对扩展方向做预条件 |
| `2` | 为初始化计算矩阵–向量乘积 |
| `3` | 为迭代更新计算矩阵–向量乘积 |

`initdvd` 和 `dvdrvr` 用列号 $j$ 指定输入及输出。`dvdson` 根据各自在工作区的起点 $p_0$，将列号转换为一基地址 $p_0+(j-1)N$。`request_slots` 的第 1 至 3 项依次给出批次列数、输入地址和输出地址。预条件时第 4 至 7 项另给出 Ritz 值、整数工作区中的待处理本征对索引、残差范数和对角预条件位移的起点。

返回时，`real_workspace` 在约束列之后连续存放这 $m$ 个谱端近似本征对的数据。数组先放 $Nm$ 个全空间向量元素，再依次放 $m$ 个本征值、$m$ 个本征值绝对变化量和 $m$ 个残差范数。低端结果按本征值升序，高端结果按降序。收敛判定针对请求态，跨度内其余位置也保留当前近似。[控制与打包源码](davidson_solver_control.f90)

接口及符号迁移见 [Davidson 本征求解库命名迁移](../../../docs/implemented/davidson_eigensolver_migration.md)。

## 来源

原 `DVDSON` 由 Andreas Stathopoulos 和 Charlotte Froese Fischer 编写。旧源码记录了 Fischer 的后续修订，并注明 `gdvd` 包装例程由 M. Saparov 编写。[原程序论文，1994](https://doi.org/10.1016/0010-4655(94)90073-6)

对角预条件和残差修正的算法背景见 Stathopoulos、Saad 与 Fischer 的论文。[Robust preconditioning，1995](https://doi.org/10.1016/0377-0427(95)00141-7)
