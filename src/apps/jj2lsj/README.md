# jj2lsj_omp：jj–LSJ 耦合变换

`jj2lsj_omp` 将 MCDHF 或 CI 得到的原子态函数从 jj 耦合基变换到 LSJ 耦合基，
用主要的 LSJ 组分标记能级。

程序读取 `isodata`、`name_desc.parquet`、`name_header.toml` 和 `name.mix` 或 `name.cimix`，
写出可供 `rlevels`、FICAL 和跃迁程序读取的 `name.mcdhf.lsj.labels` 或 `name.ci.lsj.labels`。

## 可选输出

下表的 `source` 对应 `.mix` 的 `mcdhf` 或 `.cimix` 的 `ci`。

| 文件 | 内容 |
| --- | --- |
| `name.source.lsj.csf-list` | 按转换时的对称性块列出的 LS 耦合 CSF，包括轨道占据数、子壳层项和中间耦合项。 |
| `name.source.lsj.expansion` | 所选能级的能量及带符号的 LS 展开系数。每个向量的系数顺序对应同块的 CSF 列表。 |
| `name.source.lsj.unique-labels` | 为能级分配的唯一 LSJ 标签。 |
| `name.source.lsj.unique-labels.summary` | 唯一标签分配过程的摘要。 |
| `name.source.lsj.profile.csv` | 设置 `JJ2LSJ_PROFILE=1` 时写出的分阶段计时。 |

CSF 列表供检查变换时构造的 LS 基。Suite 求解器输入使用 jj 耦合的描述符文件对。
选用非默认设置且允许省略的组分比例为零时，程序写出 CSF 列表和展开。比例大于零时，可分别选择是否写出这两份文件。
选择唯一标签时，程序写出标签及其摘要。

启用 OpenMP 后，程序并行计算 LS 系数行，线程数由 `OMP_NUM_THREADS` 控制。
写入或回读变换文件时使用串行路径。

[能级表程序](../rlevels/README.md) · [构建与运行选项](../../../docs/implemented/design-0028-jj2lsj-openmp.md)

## 计算流程

1. 装入描述符中按对称性块排列的 jj 基组与对应混合向量。根据所选能级和组分阈值选取参与构造的 jj CSF，再按壳占据和允许的耦合枚举 LSJ 基。[基组转换控制](jj2lsj_conversion.f90)
2. 对每个允许的 LSJ–jj 对，从单壳层表取变换元素，并用 9j 重耦合递推组合各壳。对能级 $s$，逐项求 $a_{is}^{\mathrm{LS}}=\sum_j T_{ij}c_{js}^{\mathrm{jj}}$。OpenMP 路径分批计算独立的 LS 行，依原有行序提交结果、累计权重及执行截断。写入或回读 `.lsj.transform.cache` 的模式走串行计算路径。[变换与并行顺序](jj2lsj_conversion.f90) · [OpenMP 实现说明](../../../docs/implemented/design-0028-jj2lsj-openmp.md)
3. 逐能级按 $a_{is}^2$ 稳定排序，保留并列组分原有顺序。按输出选项写 LS 与 jj 组分、主导标签，以及经重复标签处理的唯一标签。[组分排序](jj2lsj_postprocessing.f90) · [标签写出](jj2lsj_conversion.f90)

## 逐文件计算与数据流程

以 Suite `72eaa7b0` 为基线，按当前工作树中的应用源码和 CMake 清单核对。

### 流程与输入

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`jj2lsj_main.f90`](jj2lsj_main.f90) | 初始化计时、物理常数和角动量所需的对数阶乘表，读取 `isodata`，调用 LSJ 转换控制器并输出计时汇总。 |
| [`read_mixing_blocks.f90`](read_mixing_blocks.f90) | 验证 `CSF--MIX` 和电子数、CSF 数、块数。混合文件的原始轨道数只用于匹配检查。按块加载根编号、能量及向量并铺入统一 CSF 下标空间，保留经描述符压缩的轨道数。 |

### 矩阵与角向系数

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`jj2lsj_conversion.f90`](jj2lsj_conversion.f90) | 逐块筛选 jj 组分并枚举与占据、宇称和 $J$ 相容的 LS CSF。递归结合单壳层系数与 9j 重耦合得到 $T_{ij}$，将 jj 向量变换为 $c_i^{LS}=\sum_j T_{ij}c_j^{jj}$。OpenMP 分批算 LS 行，按原行序累计权重和截断，随后写组分及标签。回读／写变换矩阵的模式走串行路径。 |
| [`ls_jj_coefficients_f6.f90`](ls_jj_coefficients_f6.f90) | 存放 $f^6$ 的 LS–jj 变换表。每条记录给出母项和两个 jj 子壳层项的编码以及系数的符号、分子和分母，由转换器查询。 |
| [`ls_jj_coefficients_f7.f90`](ls_jj_coefficients_f7.f90) | 存放 $f^7$ 半满壳的预计算变换表。查询按母项、子壳层占据与角动量编码选取相应系数。 |
| [`ls_jj_coefficients_low_j.f90`](ls_jj_coefficients_low_j.f90) | 按总 LSJ 项和两个相对论伙伴子壳层的量子数提供 p、d 及低占据 f 壳的预计算变换记录，供转换器匹配并还原系数。 |

### 存储与输出

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`jj2lsj_postprocessing.f90`](jj2lsj_postprocessing.f90) | 在原 LS 顺序累计 $c_i^2$ 来计算范数，再用稳定的归并排序按权重降序返回索引。权重并列时保留输入次序。并行模式可按能级独立排序。 |
| [`jj2lsj_profile.f90`](jj2lsj_profile.f90) | 根据环境设置开启性能记录，用 `system_clock` 度量各转换和写文件阶段，记录块、能级、项数或字节数并报告总耗时。 |
| [`pack_ls_subshell_label.f90`](pack_ls_subshell_label.f90) | 读取子壳层占据、父项和中间耦合，按开壳／满壳规则决定需保留的占据数与标签字段，输出用于能级标记的紧凑字符串。 |

### 支持例程

| 源文件 | 计算或数据处理步骤 |
| --- | --- |
| [`ascii_digit_value.f90`](ascii_digit_value.f90) | 比较字符与按 0–9 存放的十进制字符表，命中时返回该字符的整数下标。 |
| [`orbital_letter_to_l.f90`](orbital_letter_to_l.f90) | 在谱学轨道字母序列中定位字符并换算为零基轨道角动量 $l$，供 LS 标签与子壳层处理使用。 |
