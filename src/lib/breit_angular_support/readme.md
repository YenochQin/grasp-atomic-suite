# Breit 角系数支持库

`breit_angular_support` 为双电子 Breit 角积分提供共享缓冲区的写入例程，以及 `TALK` 和 `CXK` 的外部过程接口。

| 源文件 | 职责 |
| --- | --- |
| [breit_angular_coefficient_buffer.f90](breit_angular_coefficient_buffer.f90) | 定义 `TALK` 的显式接口。`talk_fallback` 将调用方传入的角系数和径向积分标签追加到共享缓冲区。 |
| [breit_radial_integral_coefficients.f90](breit_radial_integral_coefficients.f90) | 定义 `CXK` 的显式接口。`cxk_fallback` 被调用时以 `CXK: Error` 停止。 |
| [breit_angular_coefficient_buffer_external.f90](breit_angular_coefficient_buffer_external.f90) | 将外部 `TALK` 调用转发给 `talk_fallback`。 |
| [breit_radial_integral_coefficients_external.f90](breit_radial_integral_coefficients_external.f90) | 将外部 `CXK` 调用转发给 `cxk_fallback`。 |

## 角系数缓冲

`talk_fallback` 每次调用追加一条记录，包括系数为零的情况。计数超过当前容量时，调用 `alcbuf(2)` 将容量翻倍并保留已有记录。每条记录由一个实数系数和六个整数标签组成，标签依次是左态的两个轨道索引、右态的两个轨道索引、径向积分阶数和有符号的积分类型码。 [缓冲区扩容实现](../core/angular_coefficient_buffer.f90)

调试标志 `ibug1` 非零时，例程还向 Fortran 单元 99 写出 CSF 对、轨道主量子数及对称性标签、积分阶数、类型码和系数。两个 CSF 索引仅用于这一行调试输出。

## 链接与数值实现

两个 `_external.f90` 文件各自产生一个独立的静态库成员。串行和 MPI RCI 在链接时提供自己的 `TALK` 和 `CXK`。RCI 的 `TALK` 只记录绝对值大于 `eps` 的系数，`CXK` 在应用目录中计算 Breit 有效相互作用展开中各径向积分项的系数。 [串行 RCI 的 CXK](../../apps/rci/cxk.f90) · [MPI RCI 的 CXK](../../apps/rci_mpi/cxk.f90)
