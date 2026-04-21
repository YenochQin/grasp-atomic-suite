# `gj90` 开发蓝图

## 1. 目标

目标是在现有 `rhfs90` 基础上，把每个 ASF 的 `g_J` 计算从超精细程序中分离出来，形成一个独立的新程序，暂命名为 `gj90`。新程序只负责：

- 读取 `name.c`
- 读取 `name.w`
- 读取 `name.m` 或 `name.cm`
- 计算每个 ASF 的
  - `g_J`
  - `delta g_J`
  - `total g_J`

不再承担：

- A、B 超精细常数的完整输出
- F-dependent 超精细矩阵元输出
- Zeeman 矩阵构造

## 2. 为什么优先保留 Fortran

第一版最适合继续用 Fortran 实现，原因是：

1. `rhfs90` 的核心数值链已经完整存在于 Fortran 中。
2. `g_J` 计算依赖大量模块级共享状态和已有子程序，实现语言切换会显著增加验证成本。
3. 用 Fortran 复用原始数值核，更容易做到与现有 `rhfs90` 数值逐态一致。

语言建议排序：

1. Fortran：第一版首选
2. C++：适合长期重构，不适合第一版
3. Rust：不建议第一版直接使用
4. Python：适合做驱动和批处理，不适合数值核心

## 3. `gj90` 的最小依赖链

建议复用的主链是：

```text
gj90
  -> SETCSLA
  -> GETHFD
  -> GETMIXBLOCK
  -> FACTT
  -> COMPUTE_GJ_ONLY
       -> RINTHF
       -> RINT
       -> MATELT
       -> ONEPARTICLEJJ
       -> SETQNA
```

必须保留的输入文件作用：

- `name.c`：建立 CSF 结构与耦合树
- `name.w`：提供径向轨道波函数
- `name.m/.cm`：提供 ASF 混合系数

## 4. 建议直接复用的现有模块/子程序

建议直接沿用：

- `src/lib/lib9290/setcsla.f90`
- `src/lib/lib9290/setrwfa.f90`
- `src/lib/lib9290/lodrwf.f90`
- `src/lib/lib9290/getmixblock.f90` 的 `rhfs90` 版本
- `src/appl/rhfs90/matelt.f90`
- `src/appl/rhfs90/rinthf.f90`
- `src/lib/lib9290/rint.f90`
- `src/lib/librang90/oneparticlejj.f90`
- `src/lib/lib9290/setqna.f90`

需要新建或重写的部分：

- `gj90.f90`：新入口程序
- `compute_gj_only.f90`：只计算 `g_J` 的核心例程
- `setout_gj.f90`：专门输出 `g_J` 结果的例程

## 5. `compute_gj_only` 的设计

新例程应从 `HFSGG` 精简而来，只保留 `KT=1` 路径。

### 5.1 保留内容

- `RINTGJ(i,j) = RINTHF(i,j,1)`
- `RINTDGJ(i,j) = RINT(i,j,0)`
- `GJMELT(i,j)` 和 `DGJMELT(i,j)` 的构造
- `ONEPARTICLEJJ` 生成壳层系数
- 用 `EVEC` 将 CSF 矩阵元投影到 ASF

### 5.2 删除内容

- `KT=2` 电四极分支
- `HFC(3:5,...)`
- F-dependent hyperfine matrix
- A/B 常数输出

### 5.3 重要重构点

第一版不要保留完整 `NVEC x NVEC` `GJC/DGJC` 矩阵。若目标只是每个 ASF 的 `g_J`，建议直接只累计：

- `GJC_diag(K)`
- `DGJC_diag(K)`

即对每个 ASF `K`，只保留：

```text
GJC_diag(K)   = <ASF_K | N^(1)_gJ    | ASF_K>
DGJC_diag(K)  = <ASF_K | N^(1)_delta | ASF_K>
```

然后输出：

```text
g_J       = CVAC     * GJA1 * GJC_diag(K)
delta g_J = 0.001160 * GJA1 * DGJC_diag(K)
total g_J = g_J + delta g_J
```

这会比原始 `rhfs90` 的“先算完整状态矩阵再只取对角元”更节省内存和时间。

## 6. 并行化分析

## 6.1 最主要障碍

`ONEPARTICLEJJ -> SETQNA` 这一链条依赖全局模块状态，例如：

- `m_C`
- `stat_C`
- `orb_C`

因此，第一版不适合直接做共享内存线程并行，否则容易出现线程安全问题。

## 6.2 推荐的第一版并行方案

第一版优先使用：

- MPI
- 或多进程分块并行

原因是每个进程拥有独立的模块状态副本，更容易保证正确性。

推荐切分方式：

- 对 `IC` 主循环分块
- 每个进程处理一段 `IC` 范围
- 每个进程局部累积 `GJC_diag(:)` 与 `DGJC_diag(:)`
- 最后全局归并求和

## 6.3 可安全并行的预计算部分

在主循环之前，可先并行预计算：

- `RINTGJ(i,j)`
- `RINTDGJ(i,j)`
- `GJMELT(i,j)`
- `DGJMELT(i,j)`

这部分更适合 OpenMP 或任务并行，因为它们对 `(i,j)` 对基本独立。

## 6.4 第二阶段的线程并行目标

如果未来要做 OpenMP 线程并行，需要进一步重构：

1. 把 `SETQNA` 依赖的模块数据改成显式工作区参数
2. 让 `ONEPARTICLEJJ` 使用线程私有 scratch arrays
3. 避免在主循环中写共享状态

这属于第二阶段，不建议作为第一版目标。

## 7. 分阶段开发计划

### 阶段 A：建立串行 `gj90`

目标：

- 独立运行
- 输出每个 ASF 的 `g_J`
- 数值与 `rhfs90` 对角结果一致

步骤：

1. 新建 `gj90` 入口程序
2. 接入 `SETCSLA`
3. 接入 `GETHFD`
4. 接入 `GETMIXBLOCK`
5. 接入 `FACTT`
6. 实现 `compute_gj_only`
7. 生成简洁输出文件，例如 `name.gj`

### 阶段 B：串行验证

目标：

- 用现有 `rhfs90` 输出逐态核对

建议对照项：

- `Level`
- `J`
- `Parity`
- `g_J`
- `delta g_J`
- `total g_J`

要求：

- 数值逐态一致，误差只允许来自输出格式舍入

### 阶段 C：第一版并行

目标：

- 在不重写角动量代数核心的前提下获得加速

步骤：

1. 将 `IC` 主循环按块划分
2. 每个进程读同样输入文件
3. 每个进程处理自己的 `IC` 分块
4. 局部累积 `GJC_diag` 与 `DGJC_diag`
5. 最后归并输出

### 阶段 D：进一步优化

候选优化：

- OpenMP 并行轨道对 `(i,j)` 的预计算
- 改进 `EVEC` 对角累积的缓存访问
- 减少重复调用 `ONEPARTICLEJJ`
- 将部分中间结果块化缓存

## 8. 新程序建议输出格式

建议输出文件例如 `name.gj`，内容只保留必要信息：

```text
Level   J   Parity        g_J          delta g_J       total g_J
1       2   +        ...
2       3   -        ...
```

如果后续需要兼容现有工具，也可以保留与 `.h/.ch` 接近的格式。

## 9. 最终建议

第一版推荐方案：

- 语言：Fortran
- 数值核心：直接复用 `rhfs90` 相关 Fortran 子程序
- 目标：只保留 ASF 对角 `g_J` 计算
- 并行方案：优先 MPI/多进程，不优先 OpenMP 主循环线程并行

一句话总结：

先做一个“与 `rhfs90` 数值一致的串行 `g_J-only` 程序”，再在此基础上做进程级并行，是当前最稳妥、性价比最高的路线。
