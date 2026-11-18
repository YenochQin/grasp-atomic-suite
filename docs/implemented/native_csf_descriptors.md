# Rust 原生 CSF 描述符输入

2026-11-17。原迁移基于 suite `5806884`、rCSFs `fea061f`。
本地库构建调整基于 suite `80f6e76`：reader 源码复制自 rCSFs
`f17a862d77c45164a6d31ae1e9cfeb4efbd04cf6`，现由 suite 自己保存、编译和维护。
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

`rangular_mpi` 和新启动的 `rci_mpi` 通过标准输入询问
`Use zero-first partition from descriptor header? (y/n)`。回答 `y` 时，rank 0 读取
`[block_info].zero_order_counts`，按显式块顺序广播为 `ICCUT`。每个计数必须为正且
不超过该块总数。回答 `n` 时计算完整相互作用。缺少计数的基组可用 `rcsfs zero-first`
重新生成。RCI 在新算例的速度选项之后询问分区，续算沿用 `.res` 中的设置。

原三行字符串 CSF Parquet、旧 printed-seniority/V1/归一化描述符、Int32 逻辑列、空值列均拒绝。
Parquet 为逻辑 Int8 使用物理 INT32 编码是合法的。纯闭核有明确的一个 `[0,+1]` CSF 行，
不再从 EOF 推测一个 CSF。轨道、行及显式块顺序保持不变，包括相邻同 J/宇称块。

## 代码和所有权

| 层 | 归属/实现 |
| --- | --- |
| V2 格式 | suite `src/lib/libcsf_descriptor/src/descriptor_transport.rs`；复制的 V2 契约与 rCSFs 生产者兼容 |
| Rust I/O | 同目录 `descriptor_reader.rs`：同一 header 字节快照的哈希和 TOML 解析、footer/schema/空值检查、有界批次读取、电子数/宇称/块表预扫描 |
| C ABI 1 | 同目录 `descriptor_ffi.rs`、库内 `include/rcsfs_descriptor.h`：独立句柄、明确容量的调用方缓冲、状态码及 UTF-8 错误信息；保留原符号名 |
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

reader 的完整构建源码位于 [src/lib/libcsf_descriptor](../../src/lib/libcsf_descriptor/README.md)，
包含三份 Rust 实现、C 头文件、Cargo manifest/lockfile、测试及原 MIT 许可。
CMake 将它与其他共享库一起构建，最终链接的归档为 `<build>/lib/libcsf_descriptor.a`，
安装后位于 `<prefix>/lib/libcsf_descriptor.a`。Cargo 中间产物保存在构建目录内。
构建不读取相邻 rCSFs checkout，不接收外部 reader 静态库，也不自动同步上游源码。
原有 ABI 符号名 `rcsfs_descriptor_*` 保留，零阶计数由额外的查询函数取得。
该库没有 Python、CLI 或 CSF 生成器依赖。
验证基线为 GNU Fortran 15.3.0 和 Rust/Cargo 1.99.0，锁文件选择 Arrow/Parquet 59.3.2。
压缩支持包括 ZSTD、Snappy、GZIP、Brotli、LZ4；ZSTD 后端构建使用 C 编译器。

从 suite 根目录（不安装）：

```fish
set -gx FC /usr/bin/gfortran-15
set -gx OMPI_FC /usr/bin/gfortran-15
set -gx CARGO_BUILD_JOBS 4
./configure.sh --debug --profile f2018-check --build-dir build-native
cmake --build build-native --parallel 4
```

CMake 默认用 Cargo `--locked --offline`；依赖缓存须提前准备。
缺少缓存时可用 `-DGRASP_CSF_READER_OFFLINE=OFF` 允许 Cargo 下载锁定的第三方依赖。
`GRASP_RUST_VERSION` 指定 Rust 发行版的最低版本，默认 1.99.0。
切换工具链时使用新构建目录，并通过 `GRASP_RUSTC` 和 `GRASP_CARGO` 指定对应的编译器和 Cargo。
CTest 优先使用 `rustc` 同目录的 `rustdoc`，也可通过 `GRASP_RUSTDOC` 显式指定。
旧的 `RCSFS_DESCRIPTOR_SOURCE_DIR` / `RCSFS_DESCRIPTOR_LIBRARY` 非空配置会被明确拒绝，
迁移现有缓存时移除这两项，或选择新构建目录。
生产二进制静态包含 reader，运行时不需要 Python、
Rust 工具链或 Arrow C++；仍使用原有 Fortran、MPI、BLAS 和系统运行库。
当前构建/部署验证范围为本机 Linux，未声称跨平台分发验收。

```fish
cmake --build build-native --target csf_descriptor_build --parallel 4
ctest --test-dir build-native -R '^csf_descriptor_rust$' --output-on-failure
```

开发测试使用 Python；需要原始 CSF 的测试在测试准备阶段调用独立原生 CLI，计算程序不调用它。
配置这些夹具测试时指定原生 `rcsfs` 3.1 CLI；该选项不影响库的来源或链接：

```fish
cmake -S . -B build-native -DRCSFS_TEST_CLI=/path/to/rcsfs
ctest --test-dir build-native --output-on-failure
```

旧数据可以离线转换，之后直接运行程序：

```fish
rcsfs convert original.c rcsf_descriptors.parquet
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

### 内置 reader 构建验证（2026-11-17）

基于 suite `80f6e76` 的建库调整已完成以下检查：

- GNU Fortran 15 Check，串行 + MPI + JJ2LSJ OpenMP：全量编译成功，108/108 CTest 通过。
  新增的 `csf_descriptor_rust` 包含 4 项 Rust 格式、FFI 和生命周期测试。
- Rust fmt、Clippy 全 targets（`-D warnings`）通过；矩阵驱动 8 项测试及 Python Ruff 通过。
- 仅复制 `libcsf_descriptor/` 到隔离源码目录，在没有 rCSFs 源码/库的配置中离线构建成功。
- 新 archive 与普通库同在构建 `lib/`，链接命令使用该文件；独立安装前缀下的 `lib/`、
  C 头文件和 MIT 许可均已验证，重复构建不会无故重编译。
- 旧外部 reader 库选项明确拒绝；构建不再依赖外部 reader 路径。

本次未执行完整 24 配置矩阵，也未更新或宣称通过全局历史 Fortran 告警预算。
