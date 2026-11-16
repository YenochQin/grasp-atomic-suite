# 输入状态、CPATH 与 EPS 精度修复验证

日期：2026-11-16。基准：`e47c59df12c5ad694629bf972f63f8244e245e26` 上的工作树改动。

本轮处理[静态分析报告](fortls_fortitude_analysis_2026-11-16.md)的三类优先问题。
没有执行全仓库风格修复，没有改变 MCP 文件布局、正常数值计算顺序或启用被注释的筛选条件。
构建与回归均在本地执行；没有推送或提交集群作业。

## 1. 修复内容

### MCP 输入先检查状态，再使用数据

- 四个 `setmcp.f90`：OPEN/OPENFL 和每条头记录读取后立即检查；不再先使用未成功读入的块数、
  rank 信息或基组尺寸，再累计错误。已有的块数、文件标识和 rank 布局拒绝分支统一为非零失败。
- 两个磁盘 `setcof.f90`：标签/计数、稀疏结构、系数数据以及跳过记录、REWIND 都检查状态。
  只有成功读到的 `0,0` 才表示正常块结束；EOF 和短记录不能替代哨兵。
- 两个 `setmcp_mem.f90`：计数遍和加载遍均检查全部读取；失败后不进入标签循环、分配计算或使用部分读入的数据。
  本轮合计为上述八份 MCP 读取源码的 220 个 READ/REWIND 语句位置加入或落实即时检查，另检查打开文件的状态。
- 报告点名的两份 `lodcslBio.f90` 同步修复：先检查 IOSTAT，再读取 `record(1:2)`；正常 CSF 结束的 EOF 保留，
  缺少必需的子壳层/耦合记录则报错退出。

新增公共 [input_io](../../src/lib/libmod/input_io.f90) 和 MPI 包装
[input_io_mpi](../../src/lib/mpi90/input_io_mpi.f90)。失败信息包含过程/读取阶段、unit、IOSTAT 和可获得的文件名。
已有输入一致性检查以明确的 validation failure 报告，不将状态 0 描述成 READ 失败。

串行路径 `error stop 1`；MPI 已初始化时调用 `MPI_Abort(MPI_COMM_WORLD,1,...)`，不插入 barrier。
这避免了单个 rank 读坏文件后，其他 rank 等待 collective 而错误路径又等待它们的情况。
原有 `stopmpi` 的行为未改，其他调用者不受新错误处理入口影响。

这里强化的是读取状态和已有拒绝条件，**不是一个新的完整 MCP 语义验证器**；
不声称检测了所有可被成功读取但内容不合理的索引或计数。

### CPATH 分开检查分配状态与长度

[cpath.f90](../../src/lib/mpi90/cpath.f90) 的私有 `getenv_allocating` 改为：

```fortran
if(allocated(value)) then
    if(len(value) /= length) deallocate(value)
endif
```

不依赖 `.and.` 短路来保护未分配的 deferred-length 字符串。
目录选择、rank 子目录创建和返回参数约定不变：`tmpdir` 是临时根目录，当前目录进入 rank 子目录，
`permdir` 只在 rank 0 定义。

### 四个 SETCOF 的 EPS 明确为 DOUBLE

四份声明均由 `real, parameter :: eps = 1.0d-10` 改为
`real(double), parameter :: eps = 1.0d-10`。

进一步核对发现，**当前四个过程中的 EPS 仅出现在声明和已注释的判断中**。
因此本轮是声明一致性修复，不是重新启用系数阈值；不能据此宣称提高了现有原子计算精度。
两个 memory 版本只发生这项声明变化，其 Release `setcof.f90.o` 与基准逐字节一致。

## 2. 先失败、后通过的证据

使用实际的基线 `setmcp_mem` 和生产分配器，合法块/合法空块通过，但缺少标签记录的文件仍被加载器接受并返回。
测试随后报告：

```text
valid: accepted
empty: accepted
ERROR STOP Corrupt MCP was accepted
```

该缺失标签情形在基线二进制上重复两次得到同样结论。修复后在使用数据前带上下文退出。
基线复现命令（在 suite 根目录）：

```fish
../graspkit-tools/.venv/bin/python test/f2018/check_mcp_io.py \
    build-priority-fixes-20261116/baseline-memory/driver
```

四个 EPS 的测试从生产文件提取实际声明并编译 kind/value 契约；基准首先报 `EPS is not DOUBLE`，
修复后四份均通过。EPS 没有当前可执行的筛选调用点，因此没有人为加入一个筛选算法来制造数值测试。

CPATH 的不短路保证由语言规则和基线 C161 诊断确认，并未声称基线在本机必然崩溃。
修复后 C161 消失；对其实际私有函数体编译测试了未分配、同长度、重新分配、空变量、缺失变量等情形，
并通过真实两 rank CPATH 路径测试。

## 3. 新增回归覆盖

新增十项 CTest，MCP/CSF 驱动每种完整构建配置合计覆盖 110 个有效或失效输入场景：

| 测试组 | 覆盖 |
| --- | --- |
| 四个生产 SETMCP 读取器 | 合法/空块、缺文件、缺全局/块尺寸/标签/通道头、缺基组尺寸/标志、错误 magic/rank/block id、缺标签/系数/哨兵、正 IOSTAT 的短记录 |
| 磁盘 SETCOF | 实际读取与收缩入口；控制角向权重，核对 DA 系数 `0.125`，确保正常路径未改变 |
| 两个生产内存加载器 | 两遍读取及合法空列表；在计数后截断/替换输入，检查加载遍的头、标签、系数和哨兵失败 |
| 两个 CSF 文本读取器 | 合法首个 CSF/正常 EOF、缺子壳层记录、缺耦合记录 |
| EPS/环境函数 | 四份真实声明的 DOUBLE kind/值；私有环境函数的分配与长度行为 |
| CPATH MPI | 两个 rank、含空格的临时根目录、实际工作目录及 rank 0 的永久目录约定 |

MPI 异常输入只损坏最后一个 rank 的文件，其余 rank 会进入 barrier；测试要求整个作业在超时前非零退出，
并包含正确的错误上下文和 unit。没有把 MPI 测试降为串行。

两遍加载故障通过**测试专用分配器**在最后一次通道 32 分配后替换文件，生产读取器本身原样编译。
四个头部读取器测试仍使用生产 `memory_man`，并核对缓存的标签、索引和 `0.25` 系数。
这一区分使第二遍读取的错误分支可控，同时不改变生产模块的可见性或添加测试开关。

## 4. 编译与验证结果

GNU Fortran 15.3.0；`FC`、`OMPI_FC` 均显式使用 `/usr/bin/gfortran-15`。
使用新构建目录，没有安装或覆盖已有生产程序目录。

| 验证 | 结果 |
| --- | --- |
| Release，普通串行/MPI/两种 OpenMP | 全量构建成功；最终 77/77 CTest 通过 |
| Check，同上 | 全量构建成功；最终 77/77 CTest 通过，含运行时/浮点异常检查 |
| 关闭 MPI 和 OpenMP | 新输入读取测试目标构建成功；五项相关 CTest 通过，公共输入错误处理无 MPI 依赖 |
| Fortitude C161/C181/C022 | 21 份相关文件检查通过，未执行自动修复或风格规则 |
| RMCDHF 来源检查 | 492 文件检查通过；受影响、未受影响文件及前一 manifest 的篡改负例均被拒绝 |

两套完整 CTest 的首次运行均为 76/77；唯一余项是新 CPATH 测试错误地在所有 rank 检查 `permdir`。
按生产接口的 root-only 约定修正测试后，各自单项重跑通过。生产 CPATH 的目录返回约定没有修改。

精确来源检查新增
[input-safety 快照](../reference/rmcdhf_input_safety_snapshot.json)，只覆盖十份修改的导入文件，
并绑定未改动的前一[变量重命名快照](../reference/rmcdhf_variable_rename_snapshot.json)。
未修改更早的来源记录，也未放宽未受影响文件的哈希检查。

### 告警与适用范围

历史告警预算仍因覆盖范围过时而未通过，本轮未更新或放宽它。
当前 Check 清单为 1,917 个不同源码、2,138 条告警；其中新增测试会重复编译已有内核并产生测试桩的未使用参数提示。
按已有 target 的文件/target/告警类型/消息聚合比较，没有新增计数，四个 EPS 的 REAL(8)→REAL(4) 转换告警消失。
该比较不把变动后的行号作为身份，也不意味着所有既有告警风险已解决。

本轮未运行生产规模 SCF、性能基准或物理收敛实验。正常路径的验证来自已有回归和针对生产读取过程的受控输入；
不把异常输入回归当作原子结构物理精度证明。纯风格项和报告中的其他未选风险保留。

## 5. 本地记录

忽略目录 `build-priority-fixes-20261116/` 保留：

- `baseline-memory/` 与 `logs/mcp-memory-before*.log`：基线接受缺失标签文件的复现；
- `logs/eps-before.log`、`logs/cpath-before.log`：声明精度与短路模式的基线信号；
- `release/`、`check/`、`serial-only/`：构建及 JUnit；
- `logs/ctest-*.log`：完整回归和定向重跑；
- `logs/focused-lint-final.log`、`warning-comparison.json`：定向 lint 与告警比较；
- `io-check-scope.json`、`final-validation.json`：修改位置及最终汇总。

这些构建、输入副本、原始输出和本机路径不纳入版本控制。
