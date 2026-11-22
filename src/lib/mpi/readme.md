# GRASP MPI 支持库

这里的 Fortran 模块为 GRASP 的 MPI 程序提供运行状态、输入广播、工作目录管理和稀疏本征求解所需的操作。CMake 目标为 `grasp_mpi`。以下公式对应当前源码，径向插值、正交化和核极化势调用的公共例程也附有源码链接。

| 功能 | 源码 |
| --- | --- |
| 初始化、进程状态、消息和归约 | [mpi_runtime.f90](mpi_runtime.f90)、[mpi_state.f90](mpi_state.f90)、[mpi_utilities.f90](mpi_utilities.f90) |
| 运行时间、主机和目录汇总 | [mpi_run_logging.f90](mpi_run_logging.f90) |
| 工作目录与输入错误 | [mpi_work_paths.f90](mpi_work_paths.f90)、[mpi_change_directory.f90](mpi_change_directory.f90)、[mpi_current_directory.f90](mpi_current_directory.f90)、[mpi_make_directory.f90](mpi_make_directory.f90)、[mpi_io_errors.f90](mpi_io_errors.f90) |
| 组态状态函数（CSF）完整基组、RCI 分块读取和 RMCDHF 输入 | [mpi_csf_basis_input.f90](mpi_csf_basis_input.f90)、[mpi_csf_block_input.f90](mpi_csf_block_input.f90)、[mpi_rmcdhf_csf_input.f90](mpi_rmcdhf_csf_input.f90)、[mpi_csf_collectives.f90](mpi_csf_collectives.f90) |
| 同位素、径向波函数和核极化输入 | [mpi_isotope_input.f90](mpi_isotope_input.f90)、[mpi_radial_wavefunction_file.f90](mpi_radial_wavefunction_file.f90)、[mpi_radial_wavefunction_reader.f90](mpi_radial_wavefunction_reader.f90)、[mpi_nuclear_polarization.f90](mpi_nuclear_polarization.f90) |
| 稀疏矩阵乘法和 Davidson 初始向量 | [mpi_sparse_matrix_product.f90](mpi_sparse_matrix_product.f90)、[mpi_initial_eigenvectors.f90](mpi_initial_eigenvectors.f90) |

## MPI 运行日志

`start_mpi_application` 和 `initialize_mpi_runtime_with_status` 在 MPI 初始化后由 rank 0 记录开始时间。
`finish_mpi_application` 在集体 `MPI_Finalize` 返回后，由 rank 0 打印结束时间、总墙钟时间和完成消息。
启停例程之间用 `integer(int64)` 传递开始时的时钟计数。

`Participating nodes:` 汇总主机数、MPI rank 总数和每台主机上的 rank 数。配置 rank 本地工作目录的程序还会报告
启动目录、rank 0 的串行 I/O 目录和并行 I/O 根目录。所有 rank 的路径值相同时显示一行，多个路径按相同值归并并标注
适用的 rank，例如 `0-2,4`。各 rank 的并行 I/O 文件位于对应根目录下的编号子目录中。

## 对称稀疏哈密顿量乘法

`spicmvmpi` 用一个实对称稀疏哈密顿量 $H\in\mathbb{R}^{N\times N}$ 作用于一批向量 $B\in\mathbb{R}^{N\times M}$，得到 $C=HB$。设进程数为 $P$，Fortran 列号 $j$ 从 1 开始，则第 $j$ 列分配给编号从 0 开始的 rank $(j-1)\bmod P$。每个 rank 按列连续存放上三角非零元。第 $j$ 列的末位置 $e_j$ 由 `stored_segment_end_offsets(j)` 给出，并存放 $H_{jj}$。该列在本地数组中的起始位置为

$$
s_j=
\begin{cases}
1, & j\le P,\\
e_{j-P}+1, & j>P.
\end{cases}
$$

每个 rank 的局部结果 $C^{(p)}$ 从零开始。对所属列 $j$ 和每个输入向量 $k$，先累加对角元。对同列每个已存的非对角元 $H_{ij}$（$i<j$），`dmerge` 同时累加两个对称贡献。

$$
\begin{aligned}
C^{(p)}_{jk} &\gets C^{(p)}_{jk}+H_{jj}B_{jk}
    &&\text{for the diagonal},\\
C^{(p)}_{jk} &\gets C^{(p)}_{jk}+H_{ij}B_{ik},\\
C^{(p)}_{ik} &\gets C^{(p)}_{ik}+H_{ij}B_{jk}
    &&\text{for each stored } i<j.
\end{aligned}
$$

`gdsummpi` 用 `MPI_Allreduce` 将 $NM$ 个元素逐项求和，每个 rank 得到完整的乘积。

$$
C_{ik}=\sum_{p=0}^{P-1} C^{(p)}_{ik}
      =\sum_{j=1}^{N} H_{ij}B_{jk}.
$$

`gisummpi` 将同样的逐元素求和用于整数。`gdrsummpi_root` 和 `gdmaxmpi_root` 用 `MPI_Reduce` 在 rank 0 接收逐元素求和或 $\max_p x^{(p)}_i$。[非对角元累加实现](../core/sparse_dot_and_scatter.f90)

## Davidson 初始本征向量

`iniestmpi` 取 $H$ 的左上主子矩阵 $A=H_{1:n,1:n}$，其中 $n=\min(N,n_{\mathrm{limit}})$。各 rank 将所负责列的元素写入长度为 $n(n+1)/2$ 的打包数组。上三角列内第 $i$ 行的一基索引为

$$
\operatorname{idx}(i,j)=\frac{j(j-1)}{2}+i,\qquad 1\le i\le j\le n.
$$

对这些打包数组求全局和后，各 rank 对相同的 $A$ 调用 `DSPEVX`，按本征值索引选择最低的 $h$ 个本征对。传给 LAPACK 的绝对容差取 `2*dlamch("S")`。

$$
A u_k=\lambda_k u_k,\qquad
k=1,\ldots,h,\qquad
\lambda_1\le\cdots\le\lambda_h.
$$

每个 $u_k$ 在末尾补零到长度 $N$，作为 Davidson 的起始向量。`basis_vectors` 先按 Fortran 列顺序保存向量，再保存本征值。设 $W$ 为该一维数组，则

$$
\begin{aligned}
W_{(k-1)N+i}
  &=\begin{cases}
      (u_k)_i, & 1\le i\le n,\\
      0,       & n<i\le N,
    \end{cases}\\
W_{Nh+k} &= \lambda_k,\qquad 1\le k\le h.
\end{aligned}
$$

## CSF 数组分批广播

`broadcast_descriptor_arrays` 先校验各 rank 的数组维度和可选数组。设每列包含 $n_o$ 个轨道槽位，$\ell$ 是可选的 `csf_batch_limit`。占据数和耦合编码每列各有 $n_o$ 个元素。若携带子壳层量子数编码，该数组每列有 $3n_o$ 个元素。用最大单列宽度 $w$ 限制每次 MPI 广播的元素数，其中 $I_{\max}=\operatorname{huge}(0)$ 为默认整型的最大值。

$$
w=
\begin{cases}
3n_o,& \text{with subshell quantum codes},\\
n_o, & \text{otherwise},
\end{cases}
\qquad
b=\min\!\left(8192,\left\lfloor\frac{I_{\max}}{w}\right\rfloor,\ell\right).
$$

未指定 $\ell$ 时，最小值只取前两项。所有 rank 核对 $b$ 为同一个正整数。设总列数为 $n_{\mathrm{CSF}}$，从第 $a$ 列起，每批末列 $z$ 和调用采用的 MPI 元素数为

$$
\begin{aligned}
z&=\min(n_{\mathrm{CSF}},a+b-1),&
n_{\mathrm{batch}}&=z-a+1,\\
c_{\mathrm{occ}}&=n_o n_{\mathrm{batch}},&
c_{\mathrm{quantum}}&=3n_o n_{\mathrm{batch}},&
c_{\mathrm{coupling}}&=n_o n_{\mathrm{batch}}.
\end{aligned}
$$

## 径向波函数的近原点指数、插值与正交化

`lodrwfmpi` 按主量子数 $n$ 和相对论角量子数 $\kappa$ 匹配广播的轨道记录。设核电荷数为 $Z$，原子单位下的光速为 $c$。`orbital_origin_exponents` 用以下公式，点核分支在根号内为负时调用 `stopmpi`。

$$
\gamma_\kappa=
\begin{cases}
|\kappa|,& \text{finite nucleus},\\[2pt]
\sqrt{\kappa^2-(Z/c)^2},& \text{point nucleus}.
\end{cases}
$$

随后，每个 rank 调用 `intrpq` 将径向大分量 $P(r)$ 和小分量 $Q(r)$ 插值到当前网格。对目标半径 $r$，记按插值程序选取的源节点为 $x_1,x_2,\ldots$。对任一分量 $f$，Aitken 插值表以 $T_{m,1}(r)=f(x_m)$ 起始，递推式为

$$
T_{m,j}(r)=
\frac{
 (x_m-r)T_{j-1,j-1}(r)-(x_{j-1}-r)T_{m,j-1}(r)
}{x_m-x_{j-1}},\qquad 2\le j\le m.
$$

第 $m$ 轮估计取 $f_m=T_{m,m}(r)$。两分量均非零时，分别计算相对改变量 $|(f_m-f_{mathrm{prev}})/f_m|$，其中 $f_{mathrm{prev}}$ 为上次保存的估计，初值为零。两者均小于 `accy` 时接受结果。插值最多使用 13 个节点。随后 `intrpq` 用 `rint` 的径向积分归一化两个分量。

$$
\begin{aligned}
\langle i,j\rangle
  &=\int_0^\infty \bigl[P_i(r)P_j(r)+Q_i(r)Q_j(r)\bigr]\,\mathrm{d}r,\\
(P_i,Q_i)&\gets \frac{(P_i,Q_i)}{\sqrt{\langle i,i\rangle}}.
\end{aligned}
$$

积分在源码中用离散求积实现，两个轨道的积分上限取它们当前记录的末网格索引的较小值。所有轨道都加载后，`orthsc` 按原有顺序逐一去掉当前轨道在先前同 $\kappa$ 轨道上的投影。每次用刚更新的当前轨道重新计算重叠。相同的投影系数也用于原点系数 $a_i$。

$$
(P_i,Q_i,a_i)
 \;\gets\;
(P_i,Q_i,a_i)-\langle i,j\rangle(P_j,Q_j,a_j),
\qquad j<i,\quad\kappa_j=\kappa_i.
$$

对发生投影的轨道，例程按更新后的范数重新缩放两个分量和原点系数，再从径向尾部逐点检查 $|P_i(r)|<0.01\,\mathrm{accy}$。满足条件的尾点将两个分量都置零。[插值](../core/radial_orbital_interpolation.f90) · [径向重叠](../core/weighted_radial_overlap.f90) · [Schmidt 正交化](../core/radial_orbital_orthogonalization.f90)

## 核极化势

`npol_init_mpi` 在 rank 0 调用 `npol_build_grid`，再广播电偶极 $E1$、电四极 $E2$ 及总核极化势。$K_1,K_2$ 为通道系数，$b_1,b_2$ 为各通道截断半径。启用的通道采用下式，关闭的通道贡献为零。

$$
V_{E1}(r)=-\frac{K_1}{r^4+b_1^4},\qquad
V_{E2}(r)=-\frac{K_2}{r^6+b_2^6},\qquad
V_{\mathrm{NP}}(r)=V_{E1}(r)+V_{E2}(r).
$$

内核取 $s=\max(r,b)$ 并将分母拆成尺度和无量纲和式。对幂次 $p=4$ 或 $6$，下式中 $K/s^p$ 通过二进制尾数和指数计算。

$$
V(r)=
-\frac{K/s^p}{(r/s)^p+(b/s)^p}.
$$

GRASP 的核势数组保存 $-rV$。设原库仑数组值为 $W_{\mathrm{C}}(r)$，每个 rank 用自身径向网格上的 $r$ 和接收到的极化势值，调用 `npol_compose` 计算

$$
W_{\mathrm{total}}(r)=W_{\mathrm{C}}(r)-rV_{\mathrm{NP}}(r).
$$

$r=0$ 或 $V_{\mathrm{NP}}(r)=0$ 时沿用 $W_{\mathrm{C}}(r)$。各 rank 同步状态后，才将合成值写入核势数组。[核极化势内核](../core/nuclear_polarization_potential.f90)
