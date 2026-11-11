# rwfnestimate 的 header.toml 输入与 TF 初猜保护

状态：已实施；2026-11-12，基于 suite `c47f4d4` 的本次改动。
适用目标为普通串行 `rwfnestimate`；独立的 `rwfnestimate_csfg` 不在本次范围内。
数值复现、实际测试状态和物理验证边界见
[故障与迁移记录](../reviews/rwfnestimate_nan_header_2026-11-12.md)。

## 使用方式

在包含 `isodata` 的计算目录中运行；以下可执行文件路径按安装位置调整：

```fish
# 原文本路径，默认文件名 rcsf.inp
rwfnestimate

# 显式指定文本文件
rwfnestimate --csf-text "configuration.c"

# 直接读取 rCSFs 导出的 header，不需要 rcsf.inp 或描述符数据行
rwfnestimate --header "configuration_header.toml"

rwfnestimate --help
```

两种输入选项互斥，不能重复指定。文件路径支持空格和 Unicode，不经 shell 拼接或
外部命令解释。无参数时保留原来的交互输入顺序；之后的默认网格、旧波函数读取、
Thomas–Fermi、屏蔽氢样、自定义 Z、正交化和 `rwfn.inp` 写出流程相同。
初猜方式仍通过原来的交互菜单选择。

## Header 契约

必须有 `[header_info]` 中的五个 `header_lines` 字符串及整数 `electron_count`。
例如中性 Ni 的一个最小高角动量轨道表：

```toml
[header_info]
header_lines = [
    "Core subshells:",
    "  1s   2s   2p-  2p   3s   3p-  3p",
    "Peel subshells:",
    "  3d-  3d   4s   6h-  6h   7i-  7i",
    "CSF(s):",
]
electron_count = 28
```

`electron_count` 是 **core + peel 的总电子数**。上例 core 有 18 个电子，
不能把 peel 的 10 个电子作为总数。程序直接用这个总数设置 `NELEC`，供 TF 势使用。
core/peel 轨道顺序保持不变，`SCREENPAR` 继续根据原来的 core 定义计算氢样屏蔽参数。

要求和错误处理：

- 标签行分别为 `Core subshells:`、`Peel subshells:`、`CSF(s):`。
  轨道行可以很长，不按旧文本解析器的定长字符缓冲截断。
- 轨道标签使用 `spdfghiklm`，主量子数必须大于轨道角动量，`s-` 无效。
  core/peel 内部及两者之间不得出现重复的 `(n, kappa)`。
- 空 core、空 peel 都允许，但整个轨道表不能空；轨道总数不得超过配置中的 `NNNW`。
- 电子数必须是正整数、可表示为 suite 整数，并满足
  `core 的满占据电子数 <= electron_count <= 全部列出轨道的容量之和`。
- 缺少字段、类型错误、重复 TOML key、非法轨道及容量越界均以非零状态退出。
  对旧 header 缺少电子数的情况明确提示用新版 rCSFs 重新生成，不从元素名、文件名或
  block 数量推测电子数。

这是原生 TOML 解析，支持注释、多行数组、字符串转义、字面字符串、inline table
等合法语法。`block_info` 和 `conversion_stats` 不参与径向初猜；不扫描 CSF 数据行。
输入入口见 [rwfnestimate_input.f90](../../src/appl/rwfnestimate90/rwfnestimate_input.f90)。

## 描述符和 header 的哈希关联

读取 header 不修改文件字节，不重排字段，也不追加电子数，因而不会破坏 rCSFs V2
描述符记录的 `source_header_sha256`。

该程序只消费 header，未打开描述符 Parquet，所以不会自行检查“这个 header 是否属于
某个给定的描述符”。需要配对验证的流程继续使用 rCSFs 的描述符校验/恢复接口。
若旧 header 缺少 `electron_count`，应使用支持总电子数元数据的 rCSFs 重新导出匹配的
header 和描述符；手工追加字段会改变原始文件的哈希。

## 构建与依赖

采用 suite 的 GNU Fortran 15/CMake 构建约定。普通版初猜程序需要打开串行目标：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
cmake -S . -B build-rwfn-header-release \
    -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
    -DGRASP_BUILD_SERIAL_APPS=ON \
    -DGRASP_FORTRAN_PROFILE=f2018-release \
    -DCMAKE_BUILD_TYPE=Release
cmake --build build-rwfn-header-release --target rwfnestimate -j4
ctest --test-dir build-rwfn-header-release -R rwfnestimate --output-on-failure
```

[toml-f 0.5.2](https://github.com/toml-f/toml-f/releases/tag/v0.5.2) 的解析源码随仓库
放在 `vendor/toml-f`，以 MIT 许可使用，静态链接到普通 `rwfnestimate`。
只增加私有的 `suite_tomlf` 构建目标；不引入第二份 suite 数值库。
版本、发行归档 SHA-256 和逐文件哈希见
[UPSTREAM.json](../../vendor/toml-f/UPSTREAM.json)，
[许可文本](../../vendor/toml-f/LICENSE-MIT) 随安装保留。
编译时不下载依赖，程序运行时不调用 Python、rCSFs 或外部 TOML 转换器。

macOS 验证使用本机 Homebrew GNU Fortran **15.3.0** 并显式设置 MPI 的 Fortran 后端；
该平台没有 Linux 的 `/usr/bin/gfortran-15` 路径。这是路径差异，不能用此前 GNU 16 的
上游试验代替本次 suite/GNU 15 的验证。

## 维护测试

- `rwfnestimate_Ni_high_l`、`rwfnestimate_Zn_high_l`：默认有限核网格上，各生成
  `n<=10, l<=6` 的 88 个 TF 轨道；检查有限性、正束缚能、归一化和首次 h/i 无节点。
- `rwfnestimate_header_input`：文本与 TOML 的 TF/氢样/自定义 Z 输出逐字节比较；
  空 core、纯 core、长轨道行、长/Unicode/带空格路径、TOML 语法及输入错误；
  每次读取均检查 header 哈希不变。
- 测试核数据为静态、零自旋的 Ni-58 和 Zn-64，保存于
  [fixtures](../../test/rwfnestimate/fixtures)。它们由上游 `rnucleus` 默认有限核参数生成，
  不需要在 suite 之外安装预处理程序才能运行回归测试。

这些检查覆盖输入契约和初猜数值有效性。Ni/Zn 的实际 MCDHF 活性空间优化、网格收敛和
最终能级精度不由本次测试证明。
