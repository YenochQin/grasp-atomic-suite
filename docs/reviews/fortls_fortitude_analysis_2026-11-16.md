# fortls 与 Fortitude 当前源码分析

日期：2026-11-16。源码基准：`e47c59df12c5ad694629bf972f63f8244e245e26`。

状态：只读静态分析及针对工具诊断的语法/小样本核查已完成。没有修改计算源码、应用自动修复、
重跑原子计算或提交作业。原始诊断不能直接解释为同等数量的缺陷，下面区分源码风险、规范建议和工具误报。

后续状态：本报告保留静态分析时的基准；三类优先问题已按用户要求修复，见
[修复与验证记录](priority_input_safety_fixes_2026-11-16.md)。纯风格项未批量处理。

## 1. 工具、范围与方法

- 用户通过 uv tool 安装的 `fortls 3.2.2`、`fortitude-lint 0.10.1`；后者实际命令名为 `fortitude`。
  本次直接使用已安装入口，没有更新工具或向工作区共享 Python 环境安装依赖。
- 扫描 **1,902** 份第一方 Fortran 文本：`src/` 的 1,879 个 `.f90`、一份
  `memory_man.f90_PROFILING`，以及 `test/` 的 22 份 `.f90/.F90`。
  `vendor/toml-f` 不计入第一方 lint 数量，仅在 fortls 解析依赖时按构建关系加载。
- Fortitude 使用 `--isolated --no-fix --no-unsafe-fixes`。分别保存默认规则/默认 100 列、
  默认规则/132 列、全部稳定规则/132 列和实验性 correctness 规则的独立报告。
  项目使用 GNU free-form Fortran，主要统计采用 **132 列、Fortran 2018**。
- fortls 按现有 `compile_commands.json` 与 CMake `DependInfo.cmake` 建立 **65 个实际目标上下文**，
  使用对应的源文件、库依赖与 `-D` 宏；另对 26 个未出现在编译命令清单中的文件补充上下文。
  共 **91 个上下文**，覆盖全部 1,902 份文本，没有工具加载失败。
  不把多个应用目录的同名模块/过程混入一个全仓库符号表。
- fortls 批量检查调用已安装版本的 `LangServer.file_init`、链接解析及 `get_diagnostics` 引擎；
  最小误报样本同时通过实际 `fortls --debug_diagnostics` CLI 核对。自动更新始终禁用。
- 扫描前后全部 1,902 份源码 SHA-256 一致。编译器核查的模块输出及所有小样本均写入临时/忽略目录。

“补充文件”不等于废弃代码：例如 `test/standard_grasp/interface_capacity.f90` 由 Python 测试动态编译，
没有独立 CMake 编译条目；本次又按该 Python 脚本的实际文件清单核对了它的依赖。

## 2. 数量概览

### Fortitude

| 配置 | 诊断数 | 说明 |
| --- | ---: | --- |
| 默认规则，默认 100 列 | 36,316 | 直接展示工具默认口径 |
| 默认规则，132 列 | **35,361** | 主要分析口径：src 34,847，test 514 |
| 全部稳定规则，132 列 | 55,837 | 含一次工具 panic，结果不完整 |
| 全部稳定规则，受影响文件局部绕过 S252 后 | 56,544 | 仅 jj2lsj_code.f90 未检查 S252；仍有语法解析覆盖限制 |
| correctness 规则及其 preview，132 列 | 18,114 | 补充结果，不与上表累加 |

主口径 35,361 条按工具分类如下。`Correctness` 也包含接口和编码规范，并非已确认的错误。
Fortitude 默认把很多规则的 severity 标为 `error`；这不等于 GNU 编译错误。

| 类别 | 数量 |
| --- | ---: |
| Style | 20,848 |
| Correctness | 12,048 |
| Modernisation | 1,376 |
| Obsolescent | 974 |
| Portability | 96 |
| Error / E001 | 19 |

主要规则：

| 规则 | 数量 | 如何理解 |
| --- | ---: | --- |
| S241 单/双引号约定 | 13,590 | 默认偏好双引号；不是字符串内容错误 |
| C121 未限定 `USE ... ONLY` | 4,900 | 值得治理依赖，但当前 GRASP 共享模块结构大量有意使用它 |
| S231 关键字周围空格 | 3,754 | 风格 |
| S061 END 未重复实体名称 | 1,848 | 风格与可读性 |
| C061 形参缺少 INTENT | 1,620 | 应按实际读写和接口核定，不能机械填写 |
| MOD021 旧关系运算符写法 | 1,295 | 现代化建议 |
| C003 未禁止隐式外部过程 | 1,057 | 约束未来错误的接口规范，不等于当前所有调用缺少接口 |
| C131 未显式声明可见性策略 | 989 | 模块 API 规范 |
| C001 隐式类型规则未禁用 | 976 | 模块/程序等作用域的维护风险 |
| C002 接口体未禁用隐式类型 | 894 | 接口体需单独核对，不能假设宿主的规则自动适用 |

重复的串行/MPI/内存版本会分别计数；这里没有声称这些数量都是独立根因。

### fortls

原始结果共 **57 条独立诊断**：5 Error、18 Warning、34 Information。

| 诊断 | 数量 | 核查结论 |
| --- | ---: | --- |
| Invalid parent | 3 | `block_*` 数组赋值被工具误认成 BLOCK，见下文 |
| `double` not found | 2 | 两份 vpint.f90 被自动误判为固定格式 |
| 宿主作用域名称遮蔽 | 17 | 7 条由上述错误嵌套产生；其余主要是形参/测试局部名，需要可读性审查，不是已确认功能缺陷 |
| 超过 132 列 | 1 | ris_cal.f90:359 总长 135，注释之前为 127；是行宽提示 |
| 外部 `mpi` 模块未在源码项目中找到 | 17 | MPI 由编译工具链提供，不能据此判定 MPI 未安装或链接失败 |
| 其他模块未找到 | 17 | 14 条来自未纳入当前 CMake 的旧 dspevx.f90；3 条动态测试依赖补齐后消失 |

## 3. 优先处理的源码风险

### 3.1 MCP 读取后未及时检查 IOSTAT

主口径 C181 共 **107 条候选诊断**，不能全部视为缺陷。源码复核确认以下路径值得优先处理：

- [rmcdhf/setcof.f90](../../src/appl/rmcdhf/setcof.f90)：284 行读取 `lab,ncontr` 后直接检查/解码它们，
  424 行也有类似逻辑；未先确认 `ios == 0`。
- [rmcdhf_mem/setmcp_mem.f90](../../src/appl/rmcdhf_mem/setmcp_mem.f90)：77、85 行读取标签和计数，
  循环条件依赖读入值；后续各通道数组读取也有同类状态未检查点。memory MPI 版应一起核对。
- [rbiotransform/lodcslBio.f90](../../src/appl/rbiotransform/lodcslBio.f90)：190 行读取后，
  在判断 `ios` 前先检查 `record(1:2)`，并可能进行第二次读写覆盖状态；MPI 副本同类。

对于截断/损坏文件，失败输入后的变量不应当作可靠标签、计数或 CSF 记录使用。
应在消费数据前明确处理 EOF/读取错误，并区分正常 `0,0` 哨兵与异常文件结束。
MPI 路径还要定义一致的错误退出行为。

本轮没有注入损坏 MCP 来重现挂起或错误结果，因此结论是**已定位的错误处理风险**，
不是“已发生物理结果错误”。这些典型语句在改名前基准 `c6548d6` 中已存在。

误报示例：[rangular_input.f90](../../src/appl/rangular/rangular_input.f90):80 的 `integer_status`
在 `BLOCK` 内赋值，退出 BLOCK 后第 82 行立即检查；工具的同作用域判定没有识别这条保护路径。

### 3.2 依赖逻辑短路保护未分配字符串

C161 只有一处：[cpath.f90](../../src/lib/mpi90/cpath.f90):361：

```fortran
if(allocated(value) .and. len(value) /= length) deallocate(value)
```

`value` 是 deferred-length allocatable character。Fortran 不保证 `.and.` 短路，
不能依赖左侧为假来避免右侧查询未分配对象的长度。建议把分配状态判断与长度比较分成嵌套 IF。
这属于可移植性/未定义状态访问风险，本轮没有宣称已复现崩溃；同一语句在 `c6548d6` 中已存在。

### 3.3 四个 SETCOF 的 EPS 精度声明与字面量不一致

扩展稳定规则中，C022 有 8 条；其中四个普通 RMCDHF 版本的 `setcof.f90` 都声明：

```fortran
real, parameter :: eps = 1.0d-10
```

代表位置：[rmcdhf/setcof.f90](../../src/appl/rmcdhf/setcof.f90):52，
memory 两版为第 59 行。`D` 字面量不会把左侧默认 `real` 自动提升为双精度。
在当前 GNU Fortran 15 的独立小样本中，确认：

| 项目 | 实测 |
| --- | --- |
| `storage_size(eps)` | 32 bits |
| 双精度字面量存储大小 | 64 bits |
| eps 转回 real64 后的值 | `1.00000001335143196E-10` |
| 相对 `1.0d-10` 的差异 | 约 `1.33514E-8` |

后续核定发现，四个 EPS 的使用处均已注释，当前执行路径不使用它们。
现已改为 `real(double), parameter`，只纠正声明精度，不恢复筛选条件；不能据此断言已有结果失准或精度得到提升。

也不能把全部 C022 改为 double：[linear_algebra_interfaces.f90](../../src/lib/libmod/linear_algebra_interfaces.f90):17
的 `slamch` 对应单精度 LAPACK 返回接口，必须保持 ABI。C021 扩展结果共 97 条，包含很多
`0.0/0.5/1.0/2.0` 等可精确表示的数值，也包含 `0.001`、`0.05` 等阈值；需按用途分别核查。

## 4. 不应直接当作缺陷的诊断

- C081 三处声明初始化：`genrwf::modify`、`prsrcn::cnum`、`start::mxiter`。
  初始化确实带来隐式 SAVE。`cnum`、`mxiter` 在该过程内只读，可考虑 PARAMETER；
  `modify` 是可变状态，若以后允许重复调用，需明确每次进入时是否重置。本轮未证明现有调用链因此出错。
- C182 四处多对象分配带 STAT：OpenMP 两处、复用缓存一处、描述符头一处。
  已核查它们会立即报错停止，或逐项检查 `allocated` 后清理失败缓存；不能称为四处“未处理分配失败”。
- C191 四处不可达代码含 rangular 未实现的重启分支及旧 RETURN 后代码。
  清理可以另做，但删除保护性 STOP 不等于实现重启支持。
- preview 增加了 C202 的 794 条可变模块变量、C201 的 22 条名称遮蔽及 C153 的一个空分支提示。
  它们适合长期设计审查，不适合作为当前工程的阻断条件。

## 5. 已确认的工具限制

### fortls 3.2.2

1. [csf_descriptor_input.f90](../../src/lib/lib9290/csf_descriptor_input.f90):205、208 的
   `block_lengths(...) = ...` / `block_labels(...) = ...` 被当成 BLOCK 起点，导致第 213、267、369 行
   后续过程被判为错误嵌套，同时引出 7 条伪名称遮蔽诊断。
   已用最小数组赋值样本复现；只把临时样本中的 `block_labels` 改成 `labels`，诊断就消失。
   两个样本均通过 GNU 15 的 Fortran 2018 语法检查。生产源码没有为规避工具而改名。
2. 两份 `rci{,_mpi}/vpint.f90` 自动被识别为 fixed form，五空格缩进的 USE 因而丢失，
   继而报 `double` 不存在。在内存中按真实 free-form 格式解析后，两条 Error 都消失。
3. 外部 MPI `.mod` 依赖与动态编译测试不能完全靠普通源目录扫描解决。
   为 `interface_capacity.f90` 加入 Python 测试的实际依赖清单后，三个缺失模块提示消失，
   剩余三个是合法形参遮蔽的可读性提示。

因此，这次 **5 条 fortls Error 均已确认是解析方式/解析器问题**，不是此次变量重命名失败的证据。
fortls 不是完整编译器，也不提供内存生命周期、MPI 死锁或数值正确性的完整证明。

### Fortitude 0.10.1

1. 19 个 E001 涉及 15 份跨行 FORMAT 描述符、3 份算术表达式续行及一个带条件预处理声明的测试。
   对这 19 份文件和 fortls 报 Error 的另 3 份文件，使用现有 Check 编译命令的真实 include/宏条件，
   另加 `-fsyntax-only`，共 **22 文件、28 个编译条件检查全部通过**。
   不能把它们写成 19 处真实 Fortran 语法错误。
2. E001 会让 Fortitude 丢弃部分后续 AST 诊断及修复，故这些文件的 lint 覆盖不完整。
   编译器接受源代码不等于它们已经经过 Fortitude 全规则审查。
3. 全规则检查 `jj2lsj_code.f90` 时触发内部 panic，位置为
   `useless_return.rs:419` 的区间断言。单独启用 **S252 superfluous-else-return** 可复现；
   仅对该文件绕过 S252 后，恢复了其余 708 条规则诊断。默认规则扫描没有这个 panic。
   原始失败报告和局部补扫报告分别保留，没有修改工具安装或源码来掩盖问题。

## 6. 建议的处理顺序

1. 先审查 MCP 输入错误处理和 `cpath` 的短路依赖，为异常输入补充回归。
2. 单独核定 `setcof::eps` 及非精确实常量的 kind；修复后比较数值输出，保留外部 ABI。
3. 再分程序族收紧 IMPLICIT NONE、INTENT、USE ONLY 和模块可见性，避免全仓库一次性修改接口。
4. 风格规则单独处理。当前 13,590 条引号诊断不应挤占错误处理工作的优先级。
   暂不执行全仓库 `--fix`，也不把全部 C 类规则整体屏蔽来获得“零问题”。
5. 编辑器 fortls 应按实际程序/依赖选择 source_dirs，排除历史 build 副本；
   已知误报保留最小复现，等待工具修复或单独采用已验证的局部处理。

这轮没有修复上述源码风险，没有重新运行完整 CTest 或 SCF，也没有进行物理精度/收敛验收。

## 7. 复现与本地结果

从仓库根目录，用已安装工具复现主要 Fortitude 扫描：

```fish
mkdir -p build-static-analysis-20261116/logs
fortitude --isolated check src test src/lib/libmod/memory_man.f90_PROFILING \
    --target-std f2018 --line-length 132 --no-fix --no-unsafe-fixes \
    --output-format json \
    --output-file build-static-analysis-20261116/fortitude-default-132.json
```

发现诊断时退出码为 1 属于工具正常行为，应读取报告；不能仅凭退出码认定扫描进程崩溃。
精确文件清单及所有命令保存在被 Git 忽略的 `build-static-analysis-20261116/`：

- `scan-files.txt`、`source-hashes.json`：1,902 份输入及哈希；
- `fortitude-default*.json`、`fortitude-all-stable*.json`、`fortitude-correctness-preview-132.json`：各规则口径；
- `fortitude-all-rule-explanations.json`：安装版本的规则解释；
- `fortls-contexts.json`、`fortls-results.json`、`fortls-unique-diagnostics.json`：目标上下文及完整诊断；
- `compiler-syntax-confirmation.json`：28 次编译器语法检查的参数和结果；
- `parser-probes/`、`fortls-free-form-controls.json`、`fortls-dynamic-fixture-control.json`：误报复现与对照；
- `else-return-panic.json`、`jj2lsj-all-except-S252.json`：Fortitude panic 及局部补扫；
- `eps-precision-probe.txt`、`summary.json`、`logs/`：精度小样本、汇总和日志。

新增版本控制内容仅为本报告及文档索引；原始诊断、临时模块、工具输出及本机路径不纳入提交。
