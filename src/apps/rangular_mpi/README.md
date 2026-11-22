# rangular_mpi：并行 MCDHF 角系数

`rangular_mpi` 将 Dirac–Coulomb 哈密顿量的角系数和稀疏结构计算分配给多个 MPI 进程。
程序读取 `name_desc.parquet` 和配套的 `name_header.toml`。启用零阶优先分区时，
从头文件的 `[block_info].zero_order_counts` 取得各块的零阶 CSF 数。

各 rank 在自己的工作目录写入 `mcp.mcpidx` 分片，rank 0 在启动目录写入 `name.rangular.log`。
后续 MPI MCDHF 需要使用相同的基组、进程数和临时目录布局。

[MCP 通道说明](../rangular/README.md) · [rmcdhf_mpi](../rmcdhf_mpi/README.md) · [rmcdhf_mem_mpi](../rmcdhf_mem_mpi/README.md)

[原生 CSF 输入及零阶优先分区](../../../docs/implemented/native_csf_descriptors.md)

[格式、预算与诊断选项](../../../docs/implemented/rmcdhf_mcp_io.md)

## 计算流程

1. rank 0 打开原生描述符并验证头文件，向所有 rank 广播块信息。选择零阶优先分区时使用头文件的每块零阶 CSF 数作截止，再逐块广播占据和耦合。[MPI 输入](rangular_mpi_input.f90) · [分区与通道设置](open_mcp_outputs_mpi.f90)
2. 对一块的下三角矩阵，各 rank 按进程数为步长遍历自己负责的 CSF 列，计算 T 和 $V^{(k)}$，确定本地矩阵元位置后按积分标签排序。[分布式角系数](compute_mcp_angular_coefficients_mpi.f90) · [排序](sort_mcp_coefficients.f90)
3. 各 rank 写本地带目录的 `mcp.mcpidx`，其中包含结构、积分组和进程布局。后续 MPI MCDHF 使用同一 rank 布局读取。[通道头与布局](connect_mcp_files.f90) · [稀疏记录](write_mcp_sparsity.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`configure_debug_output.f90`](configure_debug_output.f90) | 清零角系数诊断标志。启用非默认输入时读取调试选项并按需打开调试输出单元。 |
| [`configure_debug_output_mpi.f90`](configure_debug_output_mpi.f90) | 由 rank 0 读取调试设置并广播角向诊断标志，各 rank 随后使用一致的设置。 |
| [`rangular_mpi_input.f90`](rangular_mpi_input.f90) | 解析描述符路径及计时前缀，在 rank 0 打开并核查头文件。先同步输入状态，再广播块标签、计数及各块完整的 CSF 数组。各 rank 分别记录计时。 |
| [`rangular_mpi_main.f90`](rangular_mpi_main.f90) | 建立 MPI 和各 rank 工作目录，集体读取基组头和分区设置。每块分配并广播 CSF 数组，再生成各 rank 的 MCP 文件，rank 0 记录启动目录中的日志。 |
| [`read_calculation_options.f90`](read_calculation_options.f90) | 设置全相互作用或一阶限制。启用限制时从描述符 `[block_info].zero_order_counts` 取得各块截止位置并返回验证状态。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`compute_mcp_angular_coefficients_mpi.f90`](compute_mcp_angular_coefficients_mpi.f90) | 以 MPI 进程数为步长分配矩阵列，计算本 rank 的一电子和二电子角系数，筛除小项并确定稀疏位置。排序本 rank 的 T/V 贡献并写入对应的 MCP 通道。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`connect_mcp_files.f90`](connect_mcp_files.f90) | 从轨道的最大 `2j` 确定 V 通道数并分配临时通道单元。 |
| [`locate_mcp_restart_position.f90`](locate_mcp_restart_position.f90) | 扫描 MCP 记录以定位已完成的 CSF 对，根据已保留记录恢复矩阵位置和各积分通道计数。正常的 `rangular` 和 `rangular_mpi` 主程序把续算标志初始化为假。 |
| [`open_mcp_outputs_mpi.f90`](open_mcp_outputs_mpi.f90) | rank 0 读入零阶分区并集体确认成功，向所有 rank 广播限制标志及各块截止位置，再初始化本 rank 的目录分片。 |
| [`open_summary_file.f90`](open_summary_file.f90) | 把路径关联到摘要输出单元，检查打开状态。摘要内容由输入汇总例程写入。 |
| [`open_temporary_mcp_channels.f90`](open_temporary_mcp_channels.f90) | 按块号和通道编号拼接临时文件路径，在计算前连接存放稀疏结构及未排序系数的单元。 |
| [`sort_mcp_coefficients.f90`](sort_mcp_coefficients.f90) | 从临时文件读取轨道标签、系数、CSF 列和矩阵位置，对标签执行堆排序且同步移动关联字段。将相同积分标签归组后写出贡献数与记录。 |
| [`write_input_summary.f90`](write_input_summary.f90) | 用日期时间和块、CSF、轨道计数生成摘要，记录角积分的相互作用限制与相关设置。 |
| [`write_mcp_sparsity.f90`](write_mcp_sparsity.f90) | 从本块列号及稀疏位置重建每列的结束偏移，将紧凑列界和行索引写入分片，并调用稀疏统计。 |
| [`write_mcp_sparsity_summary_mpi.f90`](write_mcp_sparsity_summary_mpi.f90) | 将各 rank 的非零项计数求和到 rank 0，由其报告总计数和下三角容量。详细模式额外写调试计数及未完成提示。 |
