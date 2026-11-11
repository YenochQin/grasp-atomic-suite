# gfortran 15 的 F2018 首轮编译探测

日期：2026-11-08。基线提交：`d1260ba404833d20d709b34d9f26b05c32a3b615`。
状态：**标准选项已切换；全量构建失败，阻断位置已分类。**
本次没有修复计算实现，没有运行 CTest 或原子计算，也没有提交 Git。

此文保留首次探测的历史失败状态。后续已恢复当前配置的编译/链接，见
[分批修复验证](f2018_gfortran15_batches_2026-11-08.md)，不以新结果覆盖本轮原始诊断。

## 1. 实际修改与构建条件

在 [CMakeLists.txt](../../CMakeLists.txt) 的 GNU 分支加入 `-std=f2018`。
为隔离语言标准切换的影响，主构建保留既有 `-fno-automatic` 和
`-fallow-argument-mismatch`；没有同时引入 `-pedantic-errors`、改变存储策略或数学库。
因此这次是 F2018 模式探测，尚不是撤除所有兼容选项后的完整标准验收。

环境：GNU Fortran **15.3.0**、CMake **4.4.5**、Open MPI **5.0.11**、系统 OpenBLAS。
`FC` 与 `OMPI_FC` 都设为 `/usr/bin/gfortran-15`；包装器后端已核对。
普通程序、全部 RMCDHF 变体、MPI、CSFG、OpenMP 与测试目标均启用。

执行的主构建命令（仓库根目录）：

```sh
export FC=/usr/bin/gfortran-15
export OMPI_FC=/usr/bin/gfortran-15
cmake -S . -B build-gf15-f2018-probe \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-15 \
  -DMPI_Fortran_COMPILER=/usr/bin/mpifort \
  -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
  -DGRASP_BUILD_CSFG=ON -DGRASP_JJ2LSJ_OPENMP=ON \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DPython3_EXECUTABLE=../graspkit-tools/.venv/bin/python
cmake --build build-gf15-f2018-probe --parallel 1 -- -k
```

使用全新 build 目录、Unix Makefiles 和单任务构建，避免上一轮发现的并行 `.mod`
输出冲突；`-k` 使 Make 在失败后继续尝试独立目标。未调用 install。
实际主编译选项是 `-fno-automatic -std=f2018 -fallow-argument-mismatch -g`；
jj2lsj 按其既有逻辑去除 `-fno-automatic`，保留线程相关选项。
检查 `compile_commands.json` 确认所有条目含 `-std=f2018`。

## 2. 一次全量构建的实际结果

构建退出码 **2**。直接报告错误的文件有 **15 个（14 个生产源文件、1 个测试文件）**，
共 36 条 error 诊断，含级联错误；不能把 36 当作独立缺陷数。

| 首轮直接阻断 | 文件位置 | 根因与影响 |
| --- | --- | --- |
| 普通角动量库 | [librang90/nmtejj.f90](../../src/lib/librang90/nmtejj.f90)，50/60/74 行 | 三处算术 IF 被 F2018 删除；`rang90` 失败，阻断大部分普通计算与性质程序 |
| MPI 辅助库 | [mpi90/sys_mkdir.f90](../../src/lib/mpi90/sys_mkdir.f90)，43 行 | `SYSTEM` 不再作为可用 intrinsic 处理，IMPLICIT NONE 下报无类型；`mpiu90` 失败 |
| 四个串行 RMCDHF 变体 | `rmcdhf*` / `rmcdhf90*` 下 `setcof.f90`、`matrix.f90` 和 mem 的 `setmcp_mem.f90`，共 10 文件 | 16 处 `INT8` 调用无法按标准 intrinsic 解析；MPI 变体尚被依赖库阻挡 |
| CSFG 内存模块 | [csfg_memory_man.f90](../../src/lib/libmod_csfg/csfg_memory_man.f90)，672/855 行 | `MIN(SIZE(p),n)` 的实参整数 kind 不同；`mod_csfg` 失败 |
| CSFG 重启模块 | [symmatrix_restart_C.f90](../../src/lib/libmod_csfg/symmatrix_restart_C.f90)，14 行 | 非标准 `INTEGER*8` 声明；同样阻断 `mod_csfg` |
| QUAD 测试 | [test/lib9290_quad.f90](../../test/lib9290_quad.f90)，15 行 | `REAL*8` 定义精度探针被拒绝，后续 kind/常量/控制结构诊断为级联 |

`mod`、`9290`、`dvd90`、`mcp90` 和部分独立测试可执行文件完成构建，但完整计算程序链
没有构建成功。本轮也没有运行这些部分产物；没有将 CTest 运行失败与编译失败混为一谈。

## 3. 补充诊断：找出被公共库失败遮住的位置

只看上表会遗漏大量依赖目标。为判断后续修改范围，在不改源码的条件下，对 CMake
登记的 **3563 个生产源码编译条目、3413 个不同源文件**做逐文件 `-fsyntax-only` 诊断。
同一源码可能用于多个目标或编译条件，故条目数不等于文件数。

诊断借用上一轮同版本 gfortran 15、同份生产源码的 legacy 构建模块，逐条保留 CMake 的
include/define/存储设置，将 `.mod` 输出放到独立临时目录。这样可以绕过库依赖失败，
但**不能证明从源码干净重建、严格模块 ABI 或完整链接通过**。未启用的来源文件不在该统计内。

| 诊断配置 | 失败编译条目 | 失败的不同生产源文件 |
| --- | ---: | ---: |
| 与当前配置相同：F2018，保留参数不匹配放行 | 61 / 3563 | **59 / 3413** |
| 仅移除 `-fallow-argument-mismatch`，其他设置相同 | 145 / 3563 | **143 / 3413** |

移除放行选项新增 **84 个失败文件**，没有使原 59 个失败文件通过。新增诊断集中于
MPI choice-buffer 的旧隐式接口及 BLAS/辅助过程的类型、kind、rank 和数组序列关联。
这些是需要逐个核查的接口诊断，不等于 84 个均已证明造成错误物理结果的缺陷。

保留放行选项时，主要根因分布如下。按“文件 + 行”去重；同一行可包含多个调用，
因此行位置数也不是独立缺陷数。

| 根因 | 不同源文件数 | 报错行位置数 | 代表位置 |
| --- | ---: | ---: | --- |
| 算术 IF | 2 | 6 | 普通/CSFG 的 `nmtejj.f90` 各三处 |
| `INT8` | 28 | 44 | 各 RMCDHF 变体，CSFG 的 `setarrs_csfg/scfmpi` |
| `INTEGER*8/*4`、`REAL*8` | 9 | 20 | 普通/CSFG jj2lsj、CSFG 重启/标签/矩阵工具 |
| `SYSTEM` | 6 | 10 | MPI 路径函数、CSFG 文件工具/SCF |
| `ACCESS` | 3 | 3 | CSFG 的文件检查与工具 |
| `SIZEOF` | 1 | 20 | CSFG `report_mem.f90`；20 行合计产生 33 条此类诊断 |
| intrinsic 实参 kind 不一致 | 1 | 2 | `csfg_memory_man` 的 `MIN` |
| 先用数组维度后声明类型 | 3 | 7 | RIS 的两份 Seltzer 密度实现、CSFG `dvdson` |
| FORMAT 缺逗号 | 6 | 7 | 三个跃迁版本、RIS 串行/MPI、CSFG `dacon` |
| FORMAT 裸 `X` | 3 | 8 | CSFG accumulate/extend/split 工具 |
| 共享 DO 终止标签/非块终止语句 | 1 | 6 | CSFG `matain.f90` |

表中只列根因类别，后续“变量无类型”“FORMAT 标签未定义”等级联诊断不另列为修复任务。
生产源文件数 59/143 也不是整个迁移完成比例。

## 4. 已证实的下一层链接障碍

`sys_chdir.f90` 和 `sys_getwd.f90` 在这次主构建中可以生成目标文件，但 `nm -u` 显示：

| 目标文件 | legacy 模式引用 | F2018 模式引用 |
| --- | --- | --- |
| `sys_chdir.f90.o` | `_gfortran_chdir_i4` | `chdir_` |
| `sys_getwd.f90.o` | `_gfortran_getcwd` | `getcwd_` |

使用实际目标文件和一个仅提供错误处理 stub 的最小调用程序做链接对照：legacy
目标文件链接成功，F2018 目标文件链接失败，明确报告 `undefined reference to chdir_`
和 `getcwd_`。没有运行这个探针。因此 **只修 `sys_mkdir` 还不够**，三个平台函数要一起处理。

建议集中到标准 C 互操作的平台适配实现，保留 suite 路径和错误协议；创建目录可评估
`EXECUTE_COMMAND_LINE`，但需处理路径引用、命令启动错误和退出码。
`EXECUTE_COMMAND_LINE('cd ...')` 只改变子进程环境，不能替代当前进程的 `CHDIR`。
CSFG `ACCESS` 的读权限检查也不能仅用文件存在性检查替代。

## 5. 最需要修改的地方及建议顺序

### 第一批：解除共享库阻断，恢复测试基础

1. **普通/CSFG 两份 `nmtejj.f90`**：将六处算术 IF 改成保持等值匹配行为的标准条件分支。
   这里是整数比较、改动集中且影响面最大。先用已知表项的命中/未命中测试保护，再构建两个角动量库。
2. **MPI 的 `sys_mkdir/sys_chdir/sys_getwd`**：同时消除编译和上述链接阻断，保留 cwd、目录、长路径及失败状态语义。
   这是后续所有 MPI 程序恢复链接的前置工作，不应靠 `-fall-intrinsics` 绕过。
3. **CSFG `csfg_memory_man` 与 `symmatrix_restart_C`**：将尺寸查询直接使用目标整数 kind，
   例如 `SIZE(p,kind=LONG)`，使 `MIN` 两个实参一致。不要先用默认 kind 计算 SIZE 再转换，
   以免大数组尺寸已溢出。`INTEGER*8` 改为统一的标准 kind，核对 `.res` 记录宽度。
4. **`test/lib9290_quad.f90`**：用共享 `DOUBLE` 或标准 `iso_fortran_env` 定义精度，移除 `REAL*8` 探针；
   保持已有解析积分 `1e-12` 验收条件。这是恢复数值测试入口的低成本修改。

### 第二批：清理各计算链的标准阻断

5. **各 RMCDHF 变体的 `INT8`**：按已有 `vast_kind_param:LONG` 统一改成 `INT(value,kind=LONG)` 等标准转换。
   常规、mem、MPI、orbopt、strict/diag、CSFG 各自同步核查，保留模型与求解器差异。
   已经溢出的整数乘积不能靠最后一次转换修复，乘法和 I/O 位宽须单独检查。
6. **CSFG `SIZEOF/ACCESS/SYSTEM`、非标准 kind 和已删除 DO**：`report_mem` 中区分标量、数组、
   分配状态和字节单位，不能将 `SIZEOF` 机械替换成只给出单元素位数的 `STORAGE_SIZE`。
   普通/CSFG jj2lsj 与生成工具的 kind/循环改动分别回归，保留已有线程逐位一致契约。
7. **RIS、跃迁和 CSFG 求解器的声明/FORMAT**：先声明维度整数再使用；补 FORMAT 逗号与明确的 `1X`。
   这类改动应保持数值运算顺序和对外文本协议。优先核对 [RIS Seltzer 密度](../../src/appl/ris4/densread_seltz.f90)、
   [跃迁 csfm](../../src/appl/rtransition/csfm.f90)、[MPI Davidson](../../src/lib/libdvd90mpi/dvdson.f90)。

### 数值安全门槛：接口问题必须在恢复生产使用前处理

**撤除 `-fallow-argument-mismatch` 是独立验收条件，不以“标准语法终于能编译”替代。**

- MPI 路径逐步迁移到显式 `mpi_f08` 接口，核对 communicator/datatype/request/status、count 与缓冲区。
  本机现有 binding 对非连续数组的支持有限，需显式验证或 pack/unpack。
- BLAS/LAPACK 建立统一显式接口，核对数组形状、连续性、返回类型和 LP64/ILP64。
  [RCI iniestdm](../../src/appl/rci/iniestdm.f90) 的 DCOPY 调用，以及各 RMCDHF 的矩阵过程是重点。
- 当前仍有 **11 个 RMCDHF/CSFG 变体文件**声明 `INTEGER, EXTERNAL :: ddot`，与双精度 DDOT 返回契约不符；
  这类错误不会因加 `-std=f2018` 自动被完整发现。需要触发 CI 阻尼分支的独立测试和正确的接口。
- [lodcslmpi](../../src/lib/mpi90/lodcslmpi.f90) 等仍用 `INTEGER(BYTE)` 保存 MPI datatype handle，
  而 `mpix_bytes` 的输出是默认 INTEGER。数据元素大小和句柄大小必须分开处理。

接口诊断数量多不等于每个 MPI choice-buffer 调用都实际传错类型；应以 MPI/BLAS 契约逐项判定，
修复真缺陷与替换旧绑定分别记录。大规模索引和存储状态改造应按既有 [迁移计划](../plans/f2018_migration_gfortran15_2026-11-08.md)
分阶段执行，不能为消除编译诊断一次性全局扩宽整数、删除 SAVE 或改变 SCF 更新次序。

## 6. 本次保留的限制与证据

本次只修改构建标准和文档，`src/` 的计算实现保持不变。
下一步建议先做第 5 节第一批及 RMCDHF 的标准整数转换，再用同一配置重新全量构建；
此前没有依据宣称 F2018 迁移完成，也没有新增数值/物理结论。

完整诊断位于忽略提交的 `build-gf15-f2018-probe/audit/`：

- `configure.log`、`build.log`、`build-exit-code.txt`、`build-diagnostics.json`：一次全量构建的证据。
- `source_diagnostics.py`、`configured-results.json`、`no-argument-relaxation-results.json` 及各 summary/failures：逐条命令与补充诊断。
- `root-diagnostic-categories.json`：按根因归类的文件/行位置。
- `path-link-results.json` 和 `path-link-*.log`：CHDIR/GETCWD 实际对象的链接对照；探针未执行。
- `source-fingerprint.json`：3465 个生产 Fortran 源文件的快照指纹；这一数目包含未编入当前目标的来源文件。

逐文件诊断依赖同一源版本的旧模块，未来代码变化后不能继续盲用这些模块。上述证据与上一轮
legacy 构建的 11 个 NP/strict 运行失败属于不同问题；本轮因编译失败没有重跑那组数值测试。
