# Breit 角系数支持库

`breit_angular_support` 提供双电子 Breit 相互作用的径向积分系数计算和角系数缓冲区写入例程。

| 源文件 | 职责 |
| --- | --- |
| [breit_angular_coefficient_buffer.f90](breit_angular_coefficient_buffer.f90) | 定义外部 `collect_breit_angular_coefficient` 的显式接口和无筛选的模块存储过程。 |
| [breit_radial_integral_coefficients.f90](breit_radial_integral_coefficients.f90) | 模块过程 `calculate_breit_radial_integral_coefficients` 计算 Breit 有效相互作用展开中各径向积分的权重。 |
| [collect_breit_angular_coefficient.f90](collect_breit_angular_coefficient.f90) | 默认外部收集过程，调用 `store_breit_angular_coefficient`。 |

## 角系数缓冲

`store_breit_angular_coefficient` 每次调用追加一条记录，包括系数为零的情况。计数超过当前容量时，调用 `manage_angular_coefficient_buffer(2)` 将容量翻倍并保留已有记录。每条记录由一个实数系数和六个整数标签组成，标签依次是左态的两个轨道索引、右态的两个轨道索引、径向积分阶数和有符号的积分类型码。 [缓冲区扩容实现](../core/angular_coefficient_buffer.f90)

调试标志 `angular_coefficient_debug_flag` 非零时，例程还向角系数调试单元写出 CSF 对、轨道主量子数及对称性标签、积分阶数、类型码和系数。两个 CSF 索引仅用于这一行调试输出。

## 链接与数值实现

`collect_breit_angular_coefficient.f90` 单独编译为静态归档成员。自旋角库和 RCI 闭芯路径直接调用外部收集过程。串行和 MPI RCI 在链接时各自提供实现，只记录绝对值大于阈值的系数。两条路径都用共享的径向积分系数模块过程计算权重。
