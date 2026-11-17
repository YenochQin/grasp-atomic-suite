# Rust 原生 CSF 描述符输入

2026-11-17。实现基于 suite `5806884`、rCSFs `fea061f`；未覆盖安装目录。
rCSFs reader 实现提交为 `f17a862d77c45164a6d31ae1e9cfeb4efbd04cf6`，源码构建应使用包含该实现的版本。
生产迁移已提交为 `8cf6b1d`，rdensity 输入补全为 `80f6e76`；已完成的迁移计划由 Git 历史保存。

## 生产输入契约

全部保留的普通 CSF 消费程序直接读取 **V2 descriptor Parquet + 字节匹配的 header TOML**。
计算进程不调用 Python、PyArrow、rCSFs CLI，不生成 GCSDESC1 或临时 CSF 文本。
旧文本选项和自动回退已删除；GCSDESC1 不再是可接受的生产输入。

| 程序 | CSF 输入 |
| --- | --- |
| `rangular`、`rangular_mpi` | 默认 `rcsf_descriptors.parquet` + `rcsf_header.toml`；可用 `--csf-descriptors FILE.parquet --header FILE.toml` 显式指定 |
| `rmcdhf{,_mem,_mpi,_mem_mpi}` | `rcsf_descriptors.parquet` + `rcsf_header.toml`；原有 `isodata`、径向文件、MCP 和交互控制继续使用 |
| `rci{,_mpi}`、`gj{,_mpi}`、`rhfs{,_mpi}`、`ris4{,_mpi}`、`rdensity{,_mpi}`、JJ2LSJ | 状态名 `name` 对应 `name_descriptors.parquet` + `name_header.toml` |
| `rbiotransform{,_mpi}`、`rtransition{,_mpi,_phase}` | 初态、末态分别提供上述状态名对应的两套文件 |
| `rwfnestimate` | 默认 `rcsf_header.toml`，或 `--header FILE.toml`；只读轨道表和含 core 的 `header_info.electron_count`，不扫描 CSF payload |

直接指定 Parquet 时，省略 `--header` 会从文件名推导：`name_descriptors.parquet` 对应
`name_header.toml`；其他 `name.parquet` 对应 `name_header.toml`。MPI rangular 的相对路径以
serial I/O 目录为基准。其他应用原有状态名及 MPI CPATH 长度限制仍然有效。

`rwfnestimate` 的 header-only 检查不证明 payload 合法。旧 header 没有电子数时，需要重新生成
配套产物，不能修改 header 后继续搭配旧 Parquet：注释、空白及换行变化也会改变绑定哈希。
完整 CSF reader 可在未声明电子数的旧 header 上从实际描述符行确定电子数，并核对所有行。

原三行字符串 CSF Parquet、旧 printed-seniority/V1/归一化描述符、Int32 逻辑列、空值列均拒绝。
Parquet 为逻辑 Int8 使用物理 INT32 编码是合法的。纯闭核有明确的一个 `[0,+1]` CSF 行，
不再从 EOF 推测一个 CSF。轨道、行及显式块顺序保持不变，包括相邻同 J/宇称块。

## 代码和所有权

| 层 | 归属/实现 |
| --- | --- |
| V2 格式 | rCSFs 的 `descriptor_transport.rs`；原 `descriptor_schema` 重导出相同接口，生产者和 reader 共用源码 |
| Rust I/O | rCSFs `descriptor_reader.rs`：同一 header 字节快照的哈希和 TOML 解析、footer/schema/空值检查、有界批次读取、电子数/宇称/块表预扫描 |
| C ABI 1 | rCSFs `descriptor_ffi.rs`、`include/rcsfs_descriptor.h`：独立句柄、明确容量的调用方缓冲、状态码及 UTF-8 错误信息 |
| Fortran 适配 | [csf_rust_ffi.f90](../../src/lib/lib9290/csf_rust_ffi.f90)、[csf_descriptor_input.f90](../../src/lib/lib9290/csf_descriptor_input.f90)：`ISO_C_BINDING`、物理态/J/K 校验与补全、三数组或仅占据数组 |
| 应用状态 | `csf_basis_input`、`csf_block_input`、`rmcdhf_descriptor_input`；区分性质程序、RCI、MCP 已安装状态的生命周期 |
| 双基组 | [csf_basis_merge.f90](../../src/lib/lib9290/csf_basis_merge.f90)：保持两个轨道序列的稳定合并，直接映射数组；顺序冲突明确拒绝，不写/重读 `SLASK` |
| MPI | `csf_descriptor_collectives`、`csf_basis_mpi`、`csf_block_mpi_input`、`rmcdhf_descriptor_mpi`；root 读取，分配/输入状态先归约，按完整 CSF 列广播，错误非零 Abort |

Rust 不复制 `terms_C`，不把物理 seniority 转成打印标签再交回 Fortran。
Fortran 保留原支持态集合、`w` 限制、`2J/2K<=126` 和内核容量限制。
64 位文件计数在进入默认 INTEGER 数组前检查；它不代表数值内核支持任意大基组。

性质程序沿用剔除全空轨道的语义，按保留下来的顺序同步移动占据和子壳层量子数。
左右态使用独立输入源和明确的参考数组；合并保持各自的相对轨道次序、CSF 次序和块前缀。
应用的 mixing/径向输入仍须与其实际基组契约一致。

RMCDHF 保留 `SETMCP` 先行，核对 core、轨道数、电子数、CSF 总数、块长度、块标签及 rank 配置。
不覆盖 `csf_counts_by_block(0:...)`，不改变 DIAG/ICCUT/LFORDR，全局 CSF 数在读入期间保持总数。
长期只保存占据和每块有符号 J/宇称；完整 J/K 校验使用单行临时数组。
EAL 标准权重的 `itjpo`、`ispar` 改由块前缀定位，避免访问未分配的全基组耦合数组。

流关闭后不再读取；重开先释放旧句柄。打开失败清空 source，块读取失败关闭流并使输出数组失效。
Rust panic 不跨 C/Fortran 边界；正常格式/I/O 错误返回状态。abort/OOM 并不是可由 `catch_unwind`
恢复的错误。自有行缓冲有界；Parquet 库还管理 footer、列页和解压缓冲，不能据此宣称固定的总 RSS 上限。

## 构建和运行

推荐使用 rCSFs 的 `native-reader` 小型 crate。它直接引用三份共享格式/reader/FFI 源码，
没有 Python、CLI 或 CSF 生成器依赖。构建需要完整的 rCSFs Git 源码目录，不能只复制 `native-reader/`。
主 crate 归档不包含嵌套构建目录；预编译库模式也可采用主 crate 的无 Python `descriptor-ffi` 构建。
当前验证工具链为 GNU Fortran 15.3.0、Rust/Cargo 1.99.0，锁文件选择 Arrow/Parquet 59.3.2。
压缩支持包括 ZSTD、Snappy、GZIP、Brotli、LZ4；ZSTD 后端构建使用 C 编译器。

从 suite 根目录（不安装）：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
set -gx CARGO_BUILD_JOBS 4
./configure.sh --debug --profile f2018-check --build-dir build-native \
    -DRCSFS_DESCRIPTOR_SOURCE_DIR=/path/to/rCSFs
cmake --build build-native --parallel 4
```

源码来源必须显式指定，不自动探测兄弟仓库或隐式下载源码。CMake 默认用 Cargo `--locked --offline`；
依赖缓存须提前准备。源码构建核对 `RCSFS_RUST_VERSION`，默认 1.99.0；切换工具链应明确指定并重新验证。
也可先独立构建 reader，然后用 `-DRCSFS_DESCRIPTOR_LIBRARY=/path/to/librcsfs_descriptor.a` 链接。
这个模式不需要在 suite 构建机器上安装 Rust。生产二进制静态包含 reader，运行时不需要 Python、
Rust 工具链或 Arrow C++；仍使用原有 Fortran、MPI、BLAS 和系统运行库。
当前构建/部署验证范围为本机 Linux，未声称跨平台分发验收。

```fish
cargo build --manifest-path /path/to/rCSFs/native-reader/Cargo.toml --release --locked --offline
```

开发测试使用 Python；需要原始 CSF 的测试在测试准备阶段调用独立原生 CLI，计算程序不调用它。
配置测试时指定当前源码构建的 `rcsfs`：

```fish
cmake -S . -B build-native -DRCSFS_TEST_CLI=/path/to/rcsfs
ctest --test-dir build-native --output-on-failure
```

旧数据可以离线转换，之后直接运行程序：

```fish
rcsfs convert original.c rcsf.parquet --num-workers 1
rcsfs gen-descriptors rcsf.parquet rcsf_descriptors.parquet --header rcsf_header.toml --num-workers 1
/path/to/suite/build-native/bin/rangular
```

## 退役及产物身份

- 原 CSF 文本 reader 已从生产源码/库/应用目标移除。`test/descriptor_oracle/` 保存不安装的历史对照实现，
  包括 `csf_descriptor_verify`；它要求完成标记，不能只凭旧 reader 的退出码判断成功。
- [GCSDESC1 转换器](../../scripts/archive/prepare_csf_descriptors.py) 只作为历史离线工具和旧格式夹具保留。
  生产程序不再接受其输出。旧计时报告描述的是当时的交换格式，不能套用到当前原生 Parquet 路径。
- JJ2LSJ 等用于互操作的文本输出、普通日志，以及核数据、mixing、MCP、径向文件不属于本次删除范围。
  双基组程序之间不再需要文本 CSF 中间产物。
- header 哈希只绑定 header；MCP/mixing 仍没有完整的逐 CSF 身份摘要。现有元数据核对不能证明
  任意旧产物与当前行序一致。对照使用固定输入和程序哈希；来源不明的 MCP 应重新生成。

## 验证记录

维护测试覆盖原生 reader/C ABI、显式闭核、深层错误清理、全数组/仅占据数组精确对照，
RMCDHF 的实际 SETMCP/加载链路与标准权重，MPI 1/2/4 ranks、单列消息、worker 契约错误、root 输入错误，
以及左右态不同轨道表的多次插入、顺序冲突和无文本的完整变换/跃迁流程。

完整回归结果、基线修正、源码身份和验证边界见
[本次迁移验收记录](../reviews/native_csf_descriptor_validation_2026-11-17.md)。
