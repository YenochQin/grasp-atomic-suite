# rbiotransform_mpi：并行轨道双正交变换

`rbiotransform_mpi` 为辐射跃迁计算变换初、末态的径向轨道与混合系数，
将两套轨道变为双正交轨道。MPI 进程分担角系数计算。

程序读取 `isodata`，以及两态各自的 `name_desc.parquet`、`name_header.toml`、`name.wfn`
和 `name.mix` 或 `name.cimix`，为每个态写出 `name.bw` 和 `name.bm` 或 `name.cbm`。
若保存角向数据，`.TB` 缓存文件头会记录 MPI rank 和进程数。

[串行版本](../rbiotransform/README.md) · [rtransition_mpi](../rtransition_mpi/)

## 计算流程

1. 广播运行设置，装入两态的径向轨道及基组，按 $\kappa$ 分组计算重叠。对每组的共同方块，以 $S^{-1}$ 的 UL 因子构造双正交旋转。处理不等维的剩余轨道。[MPI 主程序](rbiotransform_mpi_main.f90) · [轨道变换](transform_orbital_sets.f90)
2. 按 rank 循环分配 CSF 行，生成子壳转移的一体角系数。缓存头包含 rank 和进程数，读取端据此校验可复用的 `.TB`。[分布式角数据](write_biorthogonal_mcp.f90) · [缓存检查](open_biorthogonal_angular_cache.f90)
3. 按两态各自的轨道变换矩阵执行混合系数补偿变换。在 CI 向量上应用占据因子和子壳跃迁展开后，由 rank 0 保存共享结果。[反变换与写出](read_biorthogonal_mcp.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

标记为「保留源码」的文件未直接列入本目录的 CMake 目标源清单。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_radial_grid.f90`](configure_radial_grid.f90) | 按核模型设置 `RNT`、`H`、`HP`、点数和精度。在非默认模式读取覆盖值，为两套轨道建立共同插值网格。 |
| [`copy_state_csf_data.f90`](copy_state_csf_data.f90) | 把共享基组的电子、CSF、轨道计数和 $n,\kappa,l,2j$ 及标签复制到指定态的专用数组，以便后续加载另一态。 |
| [`load_biorthogonal_csf_basis_mpi.f90`](load_biorthogonal_csf_basis_mpi.f90) | 根据加载阶段释放之前的占据与耦合数组，再通过 MPI 基组 reader 载入并广播当前态的描述符。 |
| [`load_merged_csf_basis.f90`](load_merged_csf_basis.f90) | 加载指定状态的原生描述符基组及块信息。按加载阶段释放旧占据／耦合数组或保留相应状态映射，供双正交角系数处理。 |
| [`rbiotransform_mpi_main.f90`](rbiotransform_mpi_main.f90) | 初始化 MPI，广播状态名和 CI 模式。建立同一网格，读两套轨道并算双正交变换，各 rank 生成局部重耦合记录，再完成态向量变换与输出。 |
| [`read_final_radial_wavefunctions.f90`](read_final_radial_wavefunctions.f90) | 同样验证并扫描末态径向文件，按末态轨道顺序保留成功匹配的记录，对其 $P,Q$ 执行末态插值。 |
| [`read_initial_radial_wavefunctions.f90`](read_initial_radial_wavefunctions.f90) | 读初态 `G92RWF`，按量子数匹配轨道，保存能量和原点参数，再把初态大小分量插值至工作网格。 |
| [`read_state_filenames.f90`](read_state_filenames.f90) | 读取初末态名，处理同名计算并询问是否继续。若只变换特定 $J$，把用户的 `2J` 转成程序的 `2J+1` 选择列表。 |
| [`read_state_mixing_coefficients.f90`](read_state_mixing_coefficients.f90) | 打开态的混合文件并按指定块取得电子、轨道和态维数、能量偏移、目标根及 CSF 向量，供反变换器使用。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`build_subshell_transformation_matrix.f90`](build_subshell_transformation_matrix.f90) | 分解轨道变换的方块，解包三角因子并形成作用在 CI 混合系数上的反变换矩阵，保存相应的工作区偏移。 |
| [`combine_upper_lower_factors.f90`](combine_upper_lower_factors.f90) | 翻转输入方阵的索引，将 UL 分解转为可复用的 LU 分解，再翻回并解包因子供两侧轨道旋转。 |
| [`factor_lower_upper.f90`](factor_lower_upper.f90) | 按列消元将输入分解为 $A=LU$，其中 $L$ 对角为 1。把严格下三角及含对角的上三角分别打包。 |
| [`invert_banded_matrix.f90`](invert_banded_matrix.f90) | 用平面旋转消去带内元素，监测最小主元和行列式因子，随后以回代计算逆矩阵并返回状态。 |
| [`invert_dense_submatrix.f90`](invert_dense_submatrix.f90) | 将矩阵工作段作为完整带宽交给带矩阵求逆例程，检查奇异状态或行列式阈值，输出方块逆矩阵。 |
| [`prepare_biorthogonal_subshells.f90`](prepare_biorthogonal_subshells.f90) | 收集两态的 $\kappa$ 值并建立对称组，生成原轨道、组内顺序与闭核／活跃分区间的映射，为方阵重叠和不等维变换准备索引。 |
| [`transform_orbital_sets.f90`](transform_orbital_sets.f90) | 各 $\kappa$ 组取共同方块重叠 $S$，把 $S^{-1}$ 作 UL 分解，得到满足 $C_i^{\mathsf T}SC_f=I$ 的轨道旋转。额外轨道以共同子空间的投影消去重叠，再按两侧变换矩阵更新 $P,Q$ 及反变换所需数据。 |
| [`write_biorthogonal_mcp.f90`](write_biorthogonal_mcp.f90) | 按 rank 分派各块的 CSF 行，生成活跃子壳之间的一体角系数并排序。写出带进程布局元数据的本地 `.TB` 记录。 |

### 径向积分与轨道

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`calculate_orbital_overlap_matrix.f90`](calculate_orbital_overlap_matrix.f90) | 在共同有效网格上对初末轨道计算 $\int(P_iP_f+Q_iQ_f)\,dr$，按两侧的轨道编号填入跨态重叠矩阵。 |
| [`calculate_subshell_overlap.f90`](calculate_subshell_overlap.f90) | 在指定 $\kappa$ 组和指定轨道范围内分别重算所需的跨态重叠，输出组内小矩阵供旋转处理。 |
| [`final_state_radial_integral.f90`](final_state_radial_integral.f90) | 对末态轨道对应用径向幂权与网格 Jacobian 求积，返回末态重叠或径向矩。 |
| [`initial_state_radial_integral.f90`](initial_state_radial_integral.f90) | 在初态两个轨道共同有效区间上对 $r^p(P_aP_b+Q_aQ_b)$ 求积，返回指定幂次的初态径向重叠。 |
| [`interpolate_final_orbital.f90`](interpolate_final_orbital.f90) | 对末态径向函数选择插值节点和误差控制，用 Aitken 递推计算两分量，更新归一化、原点系数与有效网格。 |
| [`interpolate_initial_orbital.f90`](interpolate_initial_orbital.f90) | 为初态目标径向点定位附近已知节点，以 Aitken 表同步插值 $P,Q$。监测相邻阶变化、归一化并设置初态有效尾部。 |

### 性质与后处理

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`apply_occupation_factor.f90`](apply_occupation_factor.f90) | 为指定子壳计算依赖该 CSF 占据的对角变换幂，逐能级把这个比例作用到对应的混合系数。 |
| [`apply_subshell_transfer.f90`](apply_subshell_transfer.f90) | 查找转移的轨道对和匹配的 MCP 角项，用一体激发系数与轨道变换元更新目标 CSF 向量，供反变换的非对角高阶项累加。 |
| [`read_biorthogonal_mcp.f90`](read_biorthogonal_mcp.f90) | 逐块读取原混合向量和已整理角记录，按 $\kappa$ 组应用反变换，用原根、能量和块维数写出变换后的 `.bm`／`.cbm`。MPI 的共享混合输出由 rank 0 写。 |
| [`transform_ci_coefficients.f90`](transform_ci_coefficients.f90) | 先累乘满占据惰性壳上的对角因子，再逐活跃子壳应用依赖占据的缩放及非对角转移。将逐阶转移按逆阶乘相加，构造补偿轨道旋转后的 CI 向量。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`open_biorthogonal_angular_cache.f90`](open_biorthogonal_angular_cache.f90) | 寻找 `.TB`，核查 CSF、轨道、对称组以及 MPI 版本的进程数和 rank 元数据，返回缓存能否使用的状态。 |
| [`print_matrix_symmetry.f90`](print_matrix_symmetry.f90) | 读取打包三角存储，按起止索引逐行输出矩阵，供 LU 分解诊断。 |
| [`sort_biorthogonal_mcp.f90`](sort_biorthogonal_mcp.f90) | 对轨道对打包标签执行堆排序，保持左右 CSF 与系数对应。按相同标签重组并记录每组贡献的边界。 |
| [`write_matrix.f90`](write_matrix.f90) | 按给定子矩阵范围格式化打印工作区，供重叠、分解和变换矩阵诊断使用。 |
| [`write_transformed_radial_wavefunctions.f90`](write_transformed_radial_wavefunctions.f90) | 分别按 `G92RWF` 写初、末态 `.bw`，保存变换后的大、小分量、能量、原点参数和对应的径向格点。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`check_orbital_order_stub.f90`](check_orbital_order_stub.f90)（保留源码） | 保留历史过程接口，入口直接返回。实际轨道顺序检查在 `prepare_biorthogonal_subshells` 中执行。 |
| [`combine_vectors.f90`](combine_vectors.f90) | 按两个给定比例合并等长实数组，用于累计 CI 转移展开项。 |
| [`copy_vector.f90`](copy_vector.f90) | 将指定长度的矩阵工作区或 CI 向量逐元素复制到目标区。 |
| [`fill_vector.f90`](fill_vector.f90) | 把连续工作向量设成同一输入值，为空矩阵、旋转工作区或向量累计量初始化。 |
| [`find_integer_extremum.f90`](find_integer_extremum.f90) | 根据模式在指定整数向量中返回最大或最小元素，非法模式停止并给出错误。 |
| [`multiply_submatrices.f90`](multiply_submatrices.f90) | 按行列范围乘两块实矩阵，可选用右矩阵的转置。逐个收缩被选的内积维度填入结果。 |
| [`scale_vector.f90`](scale_vector.f90) | 将整个输入向量原地乘给定比例，应用满壳归一因子或调整工作区。 |
| [`sum_integer_vector.f90`](sum_integer_vector.f90) | 对给定长度的整数连续片段求和，用于计算组内轨道总数或偏移。 |
| [`transpose_matrix.f90`](transpose_matrix.f90) | 交换输入二维数组的行、列下标，输出新矩阵供初态旋转因子转换。 |
| [`vector_inner_product.f90`](vector_inner_product.f90) | 将同长度实向量的对应分量相乘后求和，返回标量内积。 |
