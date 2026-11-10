# MCP 二进制解析与 rangular 回归对比

日期：2026-11-10。格式依据：源码基准 `e229657` 中普通串行/MPI 的
`SETMCP`、`SETMCP2/SETMCPmpi`、`SETSDA`、`SORT/SORTMEM`，以及 RMCDHF 的读入路径。

状态：已实现独立工具 [mcp_compare.py](../../scripts/mcp_compare.py)，
不修改 `rangular` 或 RMCDHF。需要 Python >= 3.11、NumPy >= 1.24；
工作区使用既有 `graspkit-tools/.venv`。不需要编译解析器，也不运行原子计算。

## 使用

下列命令从本仓库根目录执行。`reference` 和 `candidate` 是分别保存原版和修改版输出的目录，
各自保留生成这些文件时的 `rcsf.inp`、完整 `mcp.*` 或 `mcpNNN.*` 文件。

```bash
# 完整解析、检查结构，并统计各 T / V(k) 通道。
../graspkit-tools/.venv/bin/python scripts/mcp_compare.py inspect /path/to/reference

# 默认零容差，逐项比较每个 CSF 对的角系数贡献。
../graspkit-tools/.venv/bin/python scripts/mcp_compare.py compare \
  /path/to/reference /path/to/candidate --json /path/to/comparison.json

# 明确允许浮点差异；这里的阈值只是命令示例，不是物理精度验收标准。
../graspkit-tools/.venv/bin/python scripts/mcp_compare.py compare \
  /path/to/reference /path/to/candidate --atol 1e-12 --rtol 1e-10

# 导出 V(k=1) 的前 20 项为 JSON Lines，显示 block、CSF 对、轨道索引和数值。
../graspkit-tools/.venv/bin/python scripts/mcp_compare.py dump \
  /path/to/reference --unit 33 --limit 20
```

`dump` 支持 `--block`、`--label`；`--limit 0` 导出全部匹配项。
预览达到输出数量限制后仍检查所选通道的剩余记录；其他通道须用 `inspect` 检查。
原始系数导出保留重复贡献，不做合并。所有索引采用 Fortran 的 **1 起始**编号，CSF 编号在各 block 内独立。

退出码：`0` 表示有效/比较通过，`1` 表示结构或数值不同，`2` 表示输入、格式或命令错误。
`compare` 和 `inspect` 在标准输出给出 JSON；`--json` 同时保存报告，错误写入标准错误。
`--max-diffs 20` 控制报告中的差异示例数量，不截断扫描和差异计数。

比较前自动核对两个 `rcsf.inp` 的 SHA-256。若文件不同，拒绝把同一整数索引当成同一 CSF。
若某侧没有此文件，只有显式添加 `--assume-same-csfs` 才继续，并在报告中记录这个假设。
该选项不能绕过已发现的输入不一致。工具不能证明目录里的 `rcsf.inp` 确实是当初生成 MCP 的输入，
应在运行时保存输入输出配套副本。`KEYORB` 不存储在 MCP 中，默认按源码值 215 解码；
两侧必须使用相同的标签基数，改过该常量时显式传入 `--keyorb`。

## 比较语义

1. 比较电子数、轨道数、core、block 数量/标识/CSF 数量、KMAX、DIAG、ICCUT、LFORDR 等已存元数据。
   原文件的全局 ICCUT 只保存第一个 block 的值，不能据此恢复所有交互输入。
2. 解析 `mcp.30` 的 IENDC/IROW，将每个局部存储位置还原成 `(column CSF, row CSF)`。
   比较完整稀疏结构，包括没有显式角系数的对角元素。
3. 按 `(通道, block, 积分 LABEL, CSF 对)` 对齐贡献。同一对有多个贡献时，默认按值排序后逐项比较，
   保留贡献数量。输出排序和 MPI 进程数变化不会单独造成数值差异。
4. 判据为 `abs(reference-candidate) <= atol + rtol * max(abs(reference), abs(candidate))`。
   默认 `atol=rtol=0` 是数值完全一致；`+0/-0` 相等，不要求二进制字节一致。
   缺失项、额外项和结构差异均算不通过，即使缺失项的数值为零。
5. 若修改后的算法主动合并了同一 CSF 对/积分的多个贡献，可以加 `--sum-duplicates`。
   此模式按值排序后以 float64 求和，再比较每对总系数；报告仍给出原始贡献数量。
   这是不同的验收语义，不能证明每个原始贡献未变，也不能保证后续浮点装配逐位相同。

报告逐通道统计贡献数、匹配数、缺失数、超差数、最大绝对和相对误差；
差异示例定位到 block、T/V(k)、LABEL、轨道索引和两个 CSF 编号。
最大误差覆盖所有匹配项，不仅覆盖超差项；缺失项没有可定义的数值误差。
极端输入的绝对差超出 float64 表示范围时使用字符串输出，避免 JSON 的非法 Infinity。
结构差异与贡献差异分别计数，所以总差异数不等于不同 CSF 对的数量。

## 格式及边界

- `mcp.30`：全局结构头、三个公共头记录、各 block 的 MCP 头、NNONZ、IENDC/IROW。
- `mcp.31`：T 系数。`mcp.(32+k)`：V(k)，包括空通道。
- 系数 block 包含 `MCP, block, NCF, NCOEFF`；其后重复 `LABEL, NCONTR`
  和交错的 `(ICLMN, INDEX, REAL*8 coefficient)` 数组，以 `(0,0)` 结束。
- 自动识别大小端、4/8 字节记录标记和 4/8 字节 INTEGER/LOGICAL。
  支持 GNU Fortran 对大数组记录采用的连续子记录，包括在数值内部切分。
  小型头记录须未被切分；REAL 必须是 IEEE float64。
- 串行 `mcp.30…` 和 MPI `mcp000.30…` 自动发现。必须提供所有 ranks、所有通道；
  混放串行与 MPI 文件、遗留额外通道、缺 rank/文件、截断、尾随数据、计数不符、
  非法索引、重复稀疏位置、缺对角元素、无序积分标签和 NaN/Inf 都会拒绝。
- 保留标准的 MPI 循环列分配、每列对角项位于末尾的格式约定。
  不解析 `tmp.*` 未排序数据、CSFG 特有格式、RCI 其他角系数文件、REAL*16 或混合宽度 INTEGER。

## 性能与验证

大数组由 NumPy 批量读取、校验、排序和比较，不逐系数构造 Python 对象。
系数按积分标签流式处理；内存随全部稀疏结构和当前最大标签组增长，**不是恒定内存**。
多 rank 时同一通道的各 rank 文件同时打开；超大 rank 数还受到系统文件描述符限制。
如果单个积分标签本身有很大贡献列表，仍需足够内存容纳该组及排序临时数组。

2026-11-10 的软件验证：

- [回归测试](../../test/test_mcp_compare.py) 18 项通过，覆盖多 block、空通道/空 rank、
  串行与 2/4-rank 等价、8 种端序/记录标记/整数宽度组合、子记录、重排、重复贡献、
  容差、缺失/损坏数据和 CLI；GNU Fortran 15 实际生成的文件也通过读回。
- 本地已有 88 套 Zn I 输出全部解析通过，共读取 1,941 项贡献。
  仓库保存的一对 reference/candidate 小样本比较通过，但该样本的显式角系数通道为空，
  仅作为头部、结构和文件完整性检查证据。
- 合成性能测试：每侧 3,377,250 项系数，合计 117,091,056 字节输入；
  候选侧反转系数顺序。完整解析与比较约 1.07 s，含进程启动约 1.15 s，峰值 RSS 407 MiB。
  本机单次、缓存数据测量，不能推断服务器冷盘、超大实际模型或不同标签分布的性能。

测试入口：

```bash
../graspkit-tools/.venv/bin/python test/test_mcp_compare.py
# CMake 重新配置后可运行；缺 NumPy 时此可选工具测试标记为 skipped。
ctest --test-dir build -R '^mcp_compare$' --output-on-failure
```

这些证据验证解析和对比工具，不代表修改后的 `rangular` 已完成验收。
实际回归仍需用相同的 CSF/轨道顺序和交互设置分别运行原版与候选程序，再比较其完整输出。
MCP 一致说明所检查的角系数与稀疏结构一致，不建立完整 SCF/RCI 结果、活性空间收敛或物理精度。
