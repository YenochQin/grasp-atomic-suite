# rangular：MCDHF 角系数

`rangular` 按对称性块生成 Dirac–Coulomb 哈密顿量的一电子和二电子角系数及稀疏矩阵结构。
MCDHF 在优化轨道时将这些系数与径向积分结合，组装哈密顿量。

程序读取 `name_desc.parquet` 和配套的 `name_header.toml`。`mcp.30` 保存矩阵结构，
`mcp.31` 保存一电子角系数，`mcp.32` 起的通道保存各阶二电子角系数，
`name.rangular.log` 记录运行选项。后续串行 MCDHF 使用同一基组和生成的 MCP 文件。

[rmcdhf](../rmcdhf/README.md) · [rmcdhf_mem](../rmcdhf_mem/README.md) · [MPI 版本](../rangular_mpi/README.md)

[角系数计算流程](../../../docs/implemented/rangular_flow_and_state.md) · [原生 CSF 输入](../../../docs/implemented/native_csf_descriptors.md)

## 计算流程

1. 读 V3 描述符头和交互设置，保留原对称性块及块内 CSF 顺序。按块安装占据、单壳层态与累计耦合编码。[程序入口](rangular_main.f90) · [描述符与选项](rangular_input.f90)
2. 对每块遍历允许的下三角 CSF 对。一电子标量算符生成 T 系数，二电子角运算按多极阶 $k$ 生成 V 系数。可选的块内缓存按完整的 CSF 对状态匹配并还原系数及轨道标签。[角系数主循环](compute_mcp_angular_coefficients.f90) · [CSF 对复用](rangular_reuse.f90)
3. 给保留的矩阵元赋稀疏位置，将贡献按积分标签堆排序。`mcp.30` 存结构，`mcp.31` 存 T，`mcp.(32+k)` 存 $V^{(k)}$。后续 MCDHF 将每个角系数乘以相应径向积分再累加。[稀疏结构](write_mcp_sparsity.f90) · [通道排序与写出](sort_mcp_coefficients.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

标记为「保留源码」的文件未直接列入本目录的 CMake 目标源清单。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 清零角系数诊断标志。启用非默认输入时读取调试选项并按需打开调试输出单元。 |
| [`rangular_input.f90`](rangular_input.f90) | 解析显式描述符、头文件、计时和复用选项，拒绝与计算输出冲突的路径。通过原生 reader 校验和安装头信息，再按请求索引把一个块的占据、子壳层态与累计耦合装入数组。 |
| [`rangular_main.f90`](rangular_main.f90) | 解析描述符路径和交互选项，打开 MCP 输出。逐对称性块分配占据与耦合数组、读入本块并生成 T/V 通道和稀疏结构，最后关闭输入并写计时。 |
| [`read_calculation_options.f90`](read_calculation_options.f90) | 默认设置保留全部相互作用。选择一阶限制时，交互读取并记录每块的 `ICCUT`，把前缀留在变分空间，其余 CSF 按一阶设置处理。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`compute_mcp_angular_coefficients.f90`](compute_mcp_angular_coefficients.f90) | 枚举块的下三角 CSF 对，应用分区限制，生成一电子标量 T 和二电子 V 系数。截断小系数、排除另计的直接项并规范化轨道标签，给有效矩阵元编号。各通道按内存是否足够选择内存或临时文件排序。 |
| [`rangular_reuse.f90`](rangular_reuse.f90) | 把一对 CSF 的有序非空轨道、占据、完整子壳层态和累计耦合编码为查询键。完全匹配时映回当前轨道号并恢复原系数序列。未命中调用常规内核，预算不足时停止插入。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`connect_mcp_files.f90`](connect_mcp_files.f90) | 从轨道的最大 `2j` 确定 V 通道数。创建 `mcp.30` 起的二进制文件，把闭核数、块大小与标签写入结构文件，并写各通道的 MCP 标记和 rank 布局。 |
| [`locate_mcp_restart_position.f90`](locate_mcp_restart_position.f90) | 扫描 MCP 记录以定位已完成的 CSF 对，根据已保留记录恢复矩阵位置和各积分通道计数。正常的 `rangular` 和 `rangular_mpi` 主程序把续算标志初始化为假。 |
| [`open_mcp_outputs.f90`](open_mcp_outputs.f90) | 调用 `connect_mcp_files` 打开通道，读取相互作用选项，并将电子数、轨道数、CSF 数和对角／一阶设置写入 MCP 记录。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 把路径关联到摘要输出单元，检查打开状态。摘要内容由输入汇总例程写入。 |
| [`open_temporary_mcp_channels.f90`](open_temporary_mcp_channels.f90) | 按块号和通道编号拼接临时文件路径，在计算前连接存放稀疏结构及未排序系数的单元。 |
| [`rangular_profile.f90`](rangular_profile.f90) | 用 `system_clock` 对输入、内存分配、内核、排序等阶段计时，在启用时写阶段、块号、通道、记录数和存储路径类型。 |
| [`sort_mcp_coefficients.f90`](sort_mcp_coefficients.f90) | 从临时文件读取轨道标签、系数、CSF 列和矩阵位置，对标签执行堆排序且同步移动关联字段。将相同积分标签归组后写出贡献数与记录。 |
| [`sort_mcp_coefficients_in_memory.f90`](sort_mcp_coefficients_in_memory.f90) | 在已有通道数组中按打包积分标签堆排序，保持列号、稀疏位置与系数的对应关系，随后按组写出 MCP 贡献。 |
| [`try_allocate_sort_buffers.f90`](try_allocate_sort_buffers.f90) | 尝试同时分配标签、列号、矩阵位置及系数的排序缓冲。分配失败时释放部分已分配的数组，并通过状态让上层采用磁盘路径。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 用日期时间和块、CSF、轨道计数生成摘要，记录角积分的相互作用限制与相关设置。 |
| [`write_mcp_sparsity.f90`](write_mcp_sparsity.f90) | 从本块列号及稀疏位置重建每列的结束偏移，经回调写出对应系数记录，将列界与行索引写入结构通道。 |
| [`write_mcp_sparsity_summary.f90`](write_mcp_sparsity_summary.f90) | 读取本块已有的非零项计数，输出相对下三角容量 $N(N+1)/2$ 的比例信息。打开详细调试时额外重复计数并输出未完成详细诊断的提示。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`numeric_constants.f90`](numeric_constants.f90)（保留源码） | 提供角积分代码使用的实数常量，如分数和小量 `eps`。文件只声明参数。 |
