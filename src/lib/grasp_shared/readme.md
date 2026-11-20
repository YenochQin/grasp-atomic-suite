# GRASP 共享模块库

`grasp_shared` 为各 GRASP 程序提供共享状态、查找表和基础接口。每个 `.f90` 文件对应一个同名模块。库和模块的旧名对照见[命名迁移表](../../../docs/implemented/shared_library_naming.md)。

[轨道与 CSF 基组](orbital_csf_basis.f90)、[径向网格](radial_grid_state.f90)、[径向波函数](radial_wavefunctions.f90)和[本征态展开](eigenstate_expansions.f90)等模块保存由调用方填充的数据。[BLAS/LAPACK 接口](blas_lapack_interfaces.f90)声明外部例程的调用约定。下文给出库内实际执行的运算，以及读取查找表所需的编码约定。公式依据当前 suite 工作树，Git 基点为 `2d4c721`。

## 库内计算

### 每个 MPI rank 负责的 Slater 积分数

[distributed_slater_integrals.f90](distributed_slater_integrals.f90) 的 `rk_local_count(n, myid, nprocs)` 计算当前 rank 负责求值的条目数。[GENINTRK](../../apps/rci_mpi/genintrk.f90) 将 $N$ 个积分条目按 $i=1,\ldots,N$ 编号，并将第 $i$ 个条目分配给 $\operatorname{mod}(i,P)$，其中 $P$ 是 MPI 进程数。令 $N=qP+s$，$0\le s<P$，则 rank $r$ 的任务数为

$$
N_r=q+
\begin{cases}
1 & 1\le r\le s,\\
0 & \text{其余 rank}
\end{cases}
\qquad
q=\left\lfloor\frac{N}{P}\right\rfloor.
$$

余下的 $s$ 个条目从 rank 1 开始分配。比如 $N=10$、$P=4$ 时，rank 0 到 3 分别计算 $2,3,3,2$ 个条目，合计 10 个。函数要求 $N\ge0$、$P\ge1$ 且 $0\le r<P$。`GENINTRK` 在每个 rank 上仍为全部 $N$ 个条目分配缓存。

### 数组有效载荷的字节数

[array_storage_size.f90](array_storage_size.f90) 的 `storage_bytes(x)` 用 `STORAGE_SIZE(x)` 取得每个元素的位数 $b$。数组的维数为 $d$、各维长度为 $n_k$ 时，返回

$$
B=\left\lfloor\frac{b}{8}\right\rfloor\prod_{k=1}^{d}n_k.
$$

除以 8 使用整数除法。标量的空乘积为 1，零长度数组的乘积为 0。比如形状为 $3\times4$ 的 `real64` 数组占 96 个有效载荷字节。返回值使用 `int64`，统计范围为实参元素的存储。

### 数组重新分配时的数据保留范围

[array_memory_management.f90](array_memory_management.f90) 的 `ralloc` 先分配新目标。对于使用一基索引的重载，设第 $k$ 维的旧长度、新长度分别为 $s_k$ 和 $t_k$，它复制以下索引范围内的元素，再释放旧目标并重绑定指针。

$$
m_k=\min(s_k,t_k),
\qquad
A_{\mathrm{new}}(i_1,\ldots,i_d)=A_{\mathrm{old}}(i_1,\ldots,i_d)
\quad\text{对所有 } 1\le i_k\le m_k.
$$

带显式上下界 $[\ell,u]$ 的一维整数与字节重载使用旧目标的元素数 $s$ 计算复制终点。这个表达式把 $s-1$ 当作旧目标的末索引，对应旧目标从 0 开始编号。

$$
u_{\mathrm{copy}}=\min(s,u+1)-1,
\qquad
A_{\mathrm{new}}(i)=A_{\mathrm{old}}(i)
\quad\text{对所有 } \ell\le i\le u_{\mathrm{copy}}.
$$

新目标中超出复制范围的元素由调用方按需初始化。

### 编译期容量

[compile_time_capacities.f90](compile_time_capacities.f90) 定义径向数组容量 $N_{\max}=590$、轨道容量 $W_{\max}=127$ 和轨道索引打包进制 215。派生容量为

$$
N_{\mathrm{work}}=N_{\max}+10,\qquad
W_1=W_{\max}-1,\qquad W_2=W_{\max}-2.
$$

[grasp_parameters.f90](grasp_parameters.f90) 转导这些常量及其兼容别名。实际径向网格点数由 [radial_grid_state.f90](radial_grid_state.f90) 中的 `radial_point_count` 保存。

## 查找表的读取公式

### jj 子壳层项的行号与三元组

[subshell_term_table.f90](subshell_term_table.f90) 用 `itab` 保存每行项数，`jtab` 保存该行第一个三元组前的偏移，`ntab` 保存各项的 $(v,w,2J+1)$。这里 $v$ 为 seniority，$w$ 区分具有相同 seniority 和 $J$ 的项，$J$ 为子壳层总角动量。

按 [CSF 描述符读取器](../core/csf_descriptor_reader.f90) 的查表约定，设单电子角动量为 $j$、子壳层占据数为 $n$、容量为 $g=2j+1$。粒子与空穴共用一行。

$$
m=\min(n,g-n),
\qquad
r=
\begin{cases}
1 & m=0,\\
g(g-2)/8+m+1 & m>0.
\end{cases}
$$

空壳和满壳的第 1 行是 $(0,0,1)$。开壳层的第 $k$ 项取自以下位置，其中 $1\le k\le \texttt{itab}(r)$。

$$
o=\texttt{jtab}(r)+3(k-1),
\qquad
(v,w,2J+1)=
\bigl(\texttt{ntab}(o+1),\,\texttt{ntab}(o+2),\,\texttt{ntab}(o+3)\bigr).
$$

例如 $j=7/2$、$n=2$ 对应第 9 行，第一项是 `ntab(40:42) = (0, 0, 1)`，即 $v=0$、$J=0$。当前表的开壳层范围到 $2j=17$，其中 $2j>9$ 的各组只填入 $m=1$ 和 $m=2$。

### jj 项代码的低三组数字

[subshell_term_multiplicities.f90](subshell_term_multiplicities.f90) 和 [high_j_subshell_term_data.f90](high_j_subshell_term_data.f90) 保存整型项代码。调用方 [RUMTJJ](../spin_angular/jj_subshell_term_decoding.f90) 以 100 为进制读取代码 $T$ 的最低三组，得到两倍准自旋 $2Q$、两倍 seniority $2v$ 和两倍总角动量 $2J$。

$$
2Q=\left\lfloor\frac{T}{100^2}\right\rfloor\bmod100,\qquad
2v=\left\lfloor\frac{T}{100}\right\rfloor\bmod100,\qquad
2J=T\bmod100.
$$

例如 `mt(3) = 10203` 解码为 $(2Q,2v,2J)=(1,2,3)$。$2j=9$ 的独立表 `mt9` 使用外部项索引 $t=301,\ldots,306$，数组下标为 $t-300$。

### LS–jj 变换系数

[ls_jj_subshell_terms.f90](ls_jj_subshell_terms.f90) 定义 LS 子壳层项及系数记录 `ls_jj_me`。具体的系数表在 [jj2lsj 的数据模块](../../apps/jj2lsj/jj2lsj_data_1_C.f90)等文件中。[coeflsjjs](../../apps/jj2lsj/jj2lsj_code.f90) 匹配记录中的量子数，再将整数域 `factor`、`nom`、`denom` 转为实数系数。

$$
C_{\mathrm{LS},jj}=\texttt{factor}\sqrt{\frac{\texttt{nom}}{\texttt{denom}}}.
$$

例如记录 $(-1,2,3)$ 对应 $-\sqrt{2/3}$。$f^7$ 分支先分别对分子和分母开平方，再相除。

## 核极化状态的符号约定

[nuclear_polarization_state.f90](nuclear_polarization_state.f90) 的 `npol_capture` 保存 Coulomb 基线。按[初始化例程](../core/nuclear_polarization_initialization.f90)传入的约定，第 $i$ 个径向点保留 $\texttt{zz\_coul}(i)=-r_iV_{\mathrm{Coulomb}}(r_i)$。`npol_publish` 在网格长度与已捕获的基线一致时，复制调用方算好的 E1、E2 和总势。势的模型与计算公式见[核极化实现说明](../../../docs/implemented/rmcdhf_nuclear_polarization.md)。
