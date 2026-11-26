# RHFS `g_J` 计算流程报告

当前 Suite 从 `name_desc.parquet` 和 `name_header.toml` 读取 CSF，从 `name.wfn` 读取轨道，
从 `name.mix` 或 `name.cimix` 读取混合系数。下文的矩阵元说明保留上游内核简称。[输入与输出](../../README_ZH.md#输入与输出)

> 文档状态：已实施、串行代码流程说明。
> 核查日期：2026-10-05；源码基准：`8b298a6`。
> 本次核对了 RHFS 串行入口、轨道装载、HFSGG 累加和输出公式。本文不替代 gj/MPI 的完整操作说明；当前程序用法见[仓库 README](../../README_ZH.md)。
> 全部文档状态见[分类索引](../README.md)。

现行源码中的 `lande_g_radial_integrals` 对应下文的 `RINTGJ`，其积分核为
`r (P_i Q_j + Q_i P_j)`。`hyperfine_radial_integrals(1,:,:)` 对应磁偶极超精细项，
使用不同的径向幂次 `r^(-2)`。`inverse_j_norm` 对应 `GJA1`，
即 `1/sqrt(J(J+1))`。下文的其他过程名和变量名沿用报告核查基点。

## 1. 范围与结论

本报告分析 `src/apps/rhfs` 中 `g_J` 的实现流程。结论是：`rhfs` 并不是直接套用 Landé 因子闭式公式，而是先在 CSF 基底上构造与 `g_J` 相关的单电子约化矩阵元，再用混合系数 `EVEC` 投影到原子态，最后乘上 `J` 相关归一化因子与物理常数，输出 `g_J`、`delta g_J` 和 `total g_J`。

主调用链为：

```text
rhfs_main
  -> load_angular_csf_basis
  -> read_hyperfine_options
  -> read_mixing_blocks
  -> initialize_log_factorial_table
  -> calculate_hyperfine_matrix_elements
```

## 2. 程序入口与数据准备

入口程序是 `src/apps/rhfs/rhfs_main.f90`。它依次完成：

1. 读取状态名 `NAME` 和是否采用 CI 混合系数。
2. 调用 `SETSUM` 打开输出文件：
   - 非 CI: `name.h`、`name.hoffd`
   - CI: `name.ch`、`name.choffd`
3. 调用 `load_angular_csf_basis` 读取 V3 描述符文件对。
4. 调用 `read_hyperfine_options` 读取 `isodata`、建立径向网格、装载径向波函数。
5. 调用 `read_mixing_blocks` 装载本征值、本征矢和态标签。
6. 调用 `FACTT` 建立角动量代数所需的阶乘对数表。
7. 调用 `HFSGG` 执行超精细常数和 `g_J` 计算。

其中 `GETMIXBLOCK` 读入后会设置：

- `NCF = NCFTOT`
- `NVEC = NVECTOT`
- `EVEC(IC + (K-1)*NCF)` 为第 `K` 个原子态在第 `IC` 个 CSF 上的展开系数
- `IATJPO(K) = 2J+1`
- `IASPAR(K)` 为宇称标签
- `IVEC(K)` 为态序号

### 2.1 `name.wfn` 波函数文件的详细读取流程

`GETHFD` 并不自己解析波函数文件，而是在建立好统一径向网格后调用：

```text
GETHFD -> SETRWFA(TRIM(NAME)//'.wfn') -> LODRWF -> INTRPQ -> ORTHSC
```

#### 第一步：打开并校验文件头

`SETRWFA` 经共享读取器打开 `name.wfn`，校验 `RAD--WFN` 签名、轨道数、径向点数及文件长度。
[波函数格式](../../README.md#radial-wavefunction-stream)

#### 第二步：为所有目标轨道分配存储

真正的读取发生在 `LODRWF` 中。它首先根据 `CSL` 中已经定义好的轨道列表，为所有轨道分配：

- `PF(NNNP, NW)`：大分量径向函数
- `QF(NNNP, NW)`：小分量径向函数

同时初始化：

- `E(J) = -1`：表示该轨道尚未在 `.wfn` 文件中找到
- `GAMA(J)`：每个轨道的 `gamma` 参数
- `PF(:,J) = 0`，`QF(:,J) = 0`

目标轨道集合由读入的 CSF 基组确定。径向波函数文件为这些轨道提供径向函数。

#### 第三步：逐轨道扫描 `name.wfn`

`LODRWF` 之后进入循环，从 `.wfn` 文件里一条轨道一条轨道地读。每条轨道在文件中的记录结构是：

1. 轨道标签与长度信息：
   - `NPY`：主量子数
   - `NAKY`：Dirac 量子数 `kappa`
   - `EY`：该轨道能量
   - `MY`：该轨道在文件自带网格上的点数
2. 轨道函数本体：
   - `PZY`
   - `PA(1:MY)`：文件自带网格上的大分量
   - `QA(1:MY)`：文件自带网格上的小分量
3. 文件自带径向网格：
   - `RA(1:MY)`

Suite 流文件共用一列径向网格。每条轨道的大、小分量各有一列，只使用该轨道点数以内的数据。

#### 第四步：把文件中的轨道与当前计算所需轨道匹配

每读入一条 `(NPY, NAKY)` 轨道，程序就在当前所需轨道列表中查找第一个满足以下条件的轨道 `J`：

- `E(J) < 0`：说明还没有装载过
- `NPY == NP(J)`
- `NAKY == NAK(J)`

匹配成功后：

- `PZ(J) = PZY`
- `E(J) = EY`

然后调用 `INTRPQ(PA, QA, MY, RA, J, DNORM)`，把文件中的轨道函数装载到统一工作网格上。

如果 `.wfn` 文件扫描结束后，找到的轨道数 `NWIN` 小于所需轨道总数 `NW`，程序报错退出。这说明：`rhfs` 要求 `CSL` 中涉及到的每个轨道都必须在 `.wfn` 文件中存在。

#### 第五步：插值到统一径向网格

`INTRPQ` 的作用是把文件自带网格 `RA(1:MY)` 上的 `PA/QA`，插值到当前计算统一使用的网格 `R(1:N)` 上，结果写入：

- `PF(1:MF(J), J)`
- `QF(1:MF(J), J)`
- `MF(J)`：该轨道在统一网格上的有效末点

具体过程是：

1. 根据 `RA(MY)` 与当前统一网格 `R(1:N)` 的关系，找到该轨道实际覆盖到的最后一个统一网格点，记为 `MF(J)`。
2. 对每个统一网格点 `R(I)`，用 Aitken 插值从 `PA/QA` 得到：
   - `PF(I,J)`
   - `QF(I,J)`
3. 对尾部 `I > MF(J)` 的点全部置零。
4. 用 `RINT(J,J,0)` 计算范数，再把该轨道重新归一化。

这一步结束后，真正会被后续 `RINTHF` 和 `RINT` 使用的径向函数，已经是统一网格上的 `PF/QF`，不再是文件中原始的 `PA/QA/RA`。

#### 第六步：Schmidt 正交化

当所有所需轨道都装载完成后，`LODRWF` 调用 `ORTHSC` 做 Schmidt 正交化。它只对 `kappa` 相同的轨道进行正交化：

1. 计算同一 `kappa` 轨道间的重叠 `RINT(L,K,0)`
2. 从后面的轨道中减去前面轨道的投影
3. 重新归一化
4. 更新 `MF(L)`，并截去尾部非常小的数值点

因此，最终进入 `g_J` 计算的并不是 `.wfn` 文件原样存储的轨道，而是：

- 经过匹配
- 插值到统一网格
- 重新归一化
- 同 `kappa` 轨道正交化

之后得到的 `PF/QF/MF/PZ`

#### 第七步：哪些量真正进入 `g_J`

在 `rhfs` 的 `g_J` 计算里，真正被 `RINTHF` 和 `RINT` 直接使用的是：

- `PF(L,I)`
- `QF(L,I)`
- `MF(I)`
- 全局网格 `R(L)`、`RP(L)`

其中：

- `PF/QF` 是插值和归一化后的径向大、小分量
- `MF` 控制积分上限
- `R/RP` 来自 `GETHFD` 中先前建立好的统一径向网格

`PZ(J)`、`E(J)` 和 `GAMA(J)` 在波函数装载过程中也会被设置，但它们并不是 `RINTHF/RINT` 计算 `g_J` 径向积分时直接使用的核心数组。

## 3. `HFSGG` 中的预计算

`src/apps/rhfs/calculate_hyperfine_matrix_elements.f90` 是核心例程。开始阶段它为每对轨道 `(i,j)` 预计算径向积分和角向因子：

### 3.1 径向积分

- `RINTME(1,i,j) = RINTHF(i,j,-2)`  
  对应磁偶极超精细项
- `RINTGJ(i,j) = RINTHF(i,j,1)`  
  对应 `g_J` 主项
- `RINTME(2,i,j) = RINT(i,j,-3)`  
  对应电四极超精细项
- `RINTDGJ(i,j) = RINT(i,j,0)`  
  对应 `delta g_J` 修正项

其中：

- `RINTHF(i,j,k)` 在 `src/apps/rhfs/hyperfine_radial_integral.f90` 中定义，积分核为
  `r^k (P_i Q_j + Q_i P_j)`
- `RINT(i,j,k)` 在 `src/lib/core/weighted_radial_overlap.f90` 中定义，积分核为
  `r^k (P_i P_j + Q_i Q_j)`

### 3.2 角向矩阵元

`MATELT` 位于 `src/apps/rhfs/hyperfine_one_electron_element.f90`。对奇数阶 `K=1`，程序给出三组角向因子：

- `APART`
- `GJPART`
- `DGJPART`

其实现为：

```text
APART   = (κ1 + κ2)     * CLRX * OVLFAC
GJPART  = APART
DGJPART = -(κ1 + κ2 -1) * CLRX * OVLFAC
```

这里 `CLRX` 是 3-j 系数，定义在 `src/lib/core/relativistic_three_j.f90`。

## 4. CSF 基底上的矩阵元构造

当前串行入口在矩阵元计算前用 `load_angular_csf_basis` 安装组态和耦合数据。
[RHFS 入口](../../src/apps/rhfs/rhfs_main.f90#L114) · [V3 数组契约](design-0001-native-csf-descriptors.md)

对应的上游内核记法包括：

- 轨道列表 `NP/NAK/NH`
- 每个 CSF 的占据数信息 `IQA`
- 子壳层量子数信息 `JQSA`
- 耦合信息 `JCUPA`
- 每个 CSF 的总角动量和宇称

超精细矩阵元和一粒子角向系数计算使用已装入内存的 CSF 数据。

随后 `HFSGG` 双重循环遍历 CSF 对 `(IC, IR)`。对每个允许的张量耦合：

1. 用 `ONEPARTICLEJJ(KT, IPT, IC, IR, IA, IB, TSHELL)` 计算一体算符在壳层上的系数 `TSHELL`。
2. 若 `IA == IB`，遍历全部壳层，用 `TSHELL(IA)` 累加贡献。
3. 若 `IA != IB`，用双轨道项 `TSHELL(1)` 累加贡献。

形成的三个电子矩阵元累加器分别是：

- `ELEMNT`：超精细算符本身
- `ELEMNTGJ`：`g_J` 主项
- `ELEMNTDGJ`：`delta g_J` 修正项

只有当 `KT == 1` 且 `IDIFF == 0` 时，`g_J` 相关量才会累加。这表示 `rhfs` 中 `g_J` 只在同一 `J` 子空间内构造。

## 5. 从 CSF 到原子态

电子矩阵元得到后，程序再用态混合系数 `EVEC` 把它们投影到原子态基底：

- 对角块 `IR == IC` 直接用  
  `EVEC(IC+LOC1) * EVEC(IR+LOC2)`
- 非对角但同 `J` 的块 `IR != IC` 使用对称组合  
  `EVEC(IC+LOC1) * EVEC(IR+LOC2) + EVEC(IR+LOC1) * EVEC(IC+LOC2)`

最终存储到：

- `HFC(1,...)`：对角磁偶极 hfs
- `HFC(2,...)`：`J,J-1` 磁偶极 hfs
- `HFC(3..5,...)`：电四极 hfs
- `GJC(...)`：`g_J` 主项矩阵元
- `DGJC(...)`：`delta g_J` 修正矩阵元

注意：`rhfs` 最后只输出 `GJC` 和 `DGJC` 的对角元。

## 6. 最终 `g_J` 公式

在输出阶段，程序先由 `IATJPO(I)` 还原总角动量：

```text
J = (IATJPO(I) - 1) / 2
GJA1 = 1 / sqrt(J (J+1))
```

然后对每个对角态 `I == II` 计算：

```text
g_J       = CVAC     * GJA1 * GJC(NVEC*(I-1)+II)
delta g_J = 0.001160 * GJA1 * DGJC(NVEC*(I-1)+II)
total g_J = g_J + delta g_J
```

输出到 `name.h` 或 `name.ch`。

这里：

- `g_J` 是主项
- `delta g_J` 是额外修正项
- `total g_J` 是最终总值

从数值 `0.001160` 看，这一项对应电子异常磁矩修正；这是根据实现作出的判断，代码本身没有进一步注释。

## 7. 与超精细常数 A/B 的关系

`HFSGG` 同时计算 A、B 常数与 `g_J`，但它们在实现上是并行累加、最后分开换算：

- A/B 需要核自旋、核磁矩、核四极矩等量
- `g_J` 的输出公式只依赖 `CVAC`、`J` 和纯电子矩阵元 `GJC/DGJC`

因此，虽然 `rhfs` 被命名为超精细程序，`g_J` 在这里本质上是顺带输出的电子结构性质。

## 8. 两个实现细节

### 8.1 `J=0` 态不输出 `g_J`

输出前有条件：

```text
JJ == JJII .AND. JJII > 1
```

因为 `JJ = 2J+1`，所以 `J=0` 时 `JJ=1`，会被直接排除。这也避免了 `1/sqrt(J(J+1))` 的奇异性。

### 8.2 交互修改的 `C` 不直接进入最终 `g_J`

`GETHFD` 中允许用户修改速度光常数 `C`，但 `HFSGG` 在最终换算 `g_J` 时使用的是 `CVAC`，不是变量 `C`。所以从当前实现看，交互修改 `C` 不会直接改变最终写出的 `g_J`。

## 9. 总结

`rhfs` 中 `g_J` 的实现可以概括为：

1. 读取波函数和混合系数。
2. 预计算 `g_J` 主项与修正项的径向积分、角向因子。
3. 在 CSF 基底上通过 `ONEPARTICLEJJ` 构造一体矩阵元。
4. 用 `EVEC` 投影到原子态。
5. 对角态上按 `1/sqrt(J(J+1))` 做归一化。
6. 输出 `g_J`、`delta g_J` 和 `total g_J`。

换句话说，`rhfs` 的 `g_J` 是“多组态波函数期望值”意义下的结果，而不是单组态近似公式直接代入得到的数。
