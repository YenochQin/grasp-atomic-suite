# 中性原子高角动量初猜的 NaN 故障及 header 输入迁移

日期：2026-11-12。状态：故障已复现、修复及原生 TOML 输入已迁入 suite；
测试状态见末节，真实 Ni/Zn MCDHF 收敛尚未验证。
源码基准：上游参考 `grasp` 的 `488e875d`，suite 的 `c47f4d4`；
生成总电子数元数据的 rCSFs 基准为 `fea061f`。
该记录描述本次未提交的代码变更，不把基准提交本身视为包含新功能。

## 现象与试验目标

用户报告中性 Zn、Ni 扩展到 h/i 活性轨道、尤其首次引入 6h 时，经常出现：

```text
Enter the list of relativistic subshells:
SETPOT: Join set to NNNP/2 = 295
TAIL: Grid may be of insufficient extent
SETPOT: Join set to NNNP/2 = 295
TAIL: Grid may be of insufficient extent
...
TFWAVE: Unable to compute radial
```

试验需要区分网格范围不足和能量迭代失效，并验证修复搜索区间是否足以恢复有效的 TF
初猜。对照为相同输入的原程序、修复程序以及独立的屏蔽氢样初猜。
氢样初猜与 TF 使用不同模型，因此其成功仅作为控制，不作为 TF 物理解的替代证明。

设置如下：

| 项目 | 设置 |
| --- | --- |
| 元素 / 总电子数 | Ni：Z=Ne=28；Zn：Z=Ne=30 |
| 核模型 | Ni-58、Zn-64 默认有限核；质量设 0，核自旋和矩均设 0 |
| core | 1s、2s、2p-、2p、3s、3p-、3p，合计 18 电子 |
| 首配置 peel 占据 | Ni：3d-(4)3d(4)4s(2)；Zn：3d-(4)3d(6)4s(2) |
| 默认网格 | N=590，H=0.05，HP=0，RNT=2e-6/Z |
| 首轮轨道表 | 20 个轨道，包含 5g、6h、7h、7i、8i 两个精细结构分量 |
| 扩展回归 | n<=10、l<=6 的全部 88 个相对论轨道 |
| 判据 | TF 完成；能量为正且有限；波函数有限、非零、归一；首次 h/i 大分量无节点 |

原程序在 Ni/Zn 的 6h-、7i- 上复现迭代失败。更隐蔽的结果是：6h、7i 可能被标记为
`T-F` 成功，但能量为负，`p0`、平均半径和波函数包含 NaN。
若随后仅用氢样方法补齐“尚未成功”的 6h-/7i-，这些误报成功的轨道会被保留，最终文件
仍然无效。因此 `Execution complete` 或进程退出码 0 都不足以证明初猜有效。

## 数值失效链

`SOLVH` 中的 `E` 是正的束缚能大小，轨道能量相对自由电子阈值为 `-E`。
初始试探值使用核电荷 Z 的 Dirac 氢样能量，非相对论近似为 `Z^2/(2n^2)`。
而 `TFPOT` 的远程有效电荷为 `Z-Ne+1`，中性原子为 1。
高 l 轨道核区穿透较弱，因此用完整 Z 产生的试探束缚能可能远大于 TF 中的解。
Ni 的 6h 初值约 10.89 Hartree，修复后的 TF 解约 0.01396 Hartree。

向外/向内积分在连接点匹配大分量 P，再用小分量 Q 的不连续量更新能量：

\[
\Delta E = \frac{c P(r_j)[Q_{\mathrm{in}}(r_j)-Q_{\mathrm{out}}(r_j)]}
                 {\int(P^2+Q^2)\,dr}.
\]

离正确解很远时，这个局部修正可以很大。首次 h/i 的目标大分量节点数为
`n-l-1=0`，初始积分恰好返回零节点就会进入这个修正分支；节点数相同不能保证试探能量
已经接近正确解。

Ni 6h- 第一轮跟踪的实际数值（Hartree）：

| 量 | 值 |
| --- | ---: |
| E | 10.894576340494815 |
| EMIN | 0.011834327732911827 |
| EMAX | 12.964942654959476 |
| DELE | -69.411143856803307 |

旧版对负修正采用的上界更新为：

```fortran
EMAX = EMAX * (1.0D00 - 0.2D00 * ABS(DELE/E(IORB)))
```

当 `abs(DELE/E)>5` 时，上界变成负数。本例得到 `EMAX=-3.5554154572395191`，
小于仍为正的 EMIN。下轮试探值超出这个已损坏的区间，于是原本用于保护迭代的中点回退
也失败了：`E=(EMIN+EMAX)/2=-1.7717905647533037`。

`TAIL` 的衰减常数为：

\[
\beta=\sqrt{E(2-E/c^2)}.
\]

负的 E 使根号内变成负数，产生 NaN。原上游构建未启用浮点异常 trap，程序继续执行，
只可能在退出时报告 `IEEE_INVALID_FLAG`。这里的 invalid 与微小尾部造成的正常数值
underflow 是不同情况，不能把两者都视为可以忽略的输出。

两条警告分别来自这条链上的不同位置：

1. `SETPOT` 用试探能量寻找连接点；负 E 使本例的搜索条件无法满足，退到
   `NNNP/2=295`。295 是索引，并非物理半径或轨道编号。
2. `TAIL` 搜索满足 `BETA*(R(I)-RJP)>0.25*EXPMAX` 的尾部位置。BETA 已是 NaN，
   比较不会为真，遍历到网格末端后发出通用的范围不足警告。

Ni 此网格末端约为 `440397.1 a0`。上述警告是本次负能量失效的后果；扩大网格本身
不能修复负数开平方。其他短网格、点核或自定义设置仍可能存在真正的范围不足，不能仅凭
同一句日志就把所有算例归为这一根因。本次没有修改共享 `SETPOT` 的备用索引逻辑。

## 为什么 NaN 还能被判为成功

NaN 经尾部初始化、向内积分和连接点缩放，进入 P、Q、归一化积分和连接误差。
旧版只用下列条件决定是否继续迭代：

```fortran
IF (ABS(DMSMCH/QFJPO) > EPSLON) THEN
   ! 修正能量并重试
ENDIF
! 否则归一化，保持 FAIL=.FALSE. 并返回
```

在此次 IEEE 浮点环境中，NaN 与正常数的有序大小比较为假。
所以“误差不大于阈值这个分支未触发”不能推导出“误差小于或等于阈值”。
归一化继续产生 NaN，`FAIL` 却仍是 false，`FRMTFP` 据此设置 `SET=.TRUE.` 和 `SOURCE='T-F'`。

原版 Ni 跟踪中，6h-（kappa=+5）在失效后返回节点数 2，而目标是 0，因而迭代到上限后
报错；6h（kappa=-6）返回 0，穿过节点检查后又绕过上述误差检查，误报成功。
对含 NaN 的数组计节点没有物理意义；这里的 0/2 是此次执行的观察，不是可依赖的性质。
原输出中的代表性数据为：

```text
Shell      e           p0        gamma        <r>      MTP  SRC
  6h  -0.1769D+01         NaN  0.6000D+01         NaN  590  T-F
  7i  -0.2940D+01         NaN  0.7000D+01         NaN  590  T-F
```

## 修复与仓库边界

[SOLVH](../../src/appl/rwfnestimate/solvh.f90) 保留 TF 势、节点目标、精度阈值和
迭代上限，修改搜索边界更新：

```fortran
EMAX = MAX(E(IORB), EMAX*(1.0D00 - 0.2D00*ABS(DELE/E(IORB))))
EMIN = MIN(E(IORB), EMIN*(1.0D00 + 0.2D00*ABS(DELE/E(IORB))))
```

两式分别用于负修正和非负修正。它们保留当前合法试探值作为新边界的约束，使超出区间的
大步长可以通过原有中点逻辑重试。进入尾部前拒绝无效/非正束缚能，接受结果前检查范数、
连接误差和修正的有限性；误差判断改为乘法形式，避免除以零的小分量。

修复、回归和文档归属 suite 普通 `rwfnestimate`。先前用于试验的 `grasp` 源码修改及
新测试已经撤回，原安装二进制保持原版，构建目录中的 `rwfnestimate` 也已用恢复的源码
重新编译。没有把 suite 的共享参数、库或 Git 元数据复制到上游参考目录。

同时新增 [header.toml 输入](../implemented/rwfnestimate_header_input.md)，从五行轨道
表和含 core 的总电子数安装状态，然后进入相同的数值流程。两种输入方式通过输出字节
一致性进行对照；读 header 不改变与描述符关联的文件哈希。

## 验证证据及限制

上游阶段采用 GNU 16.2.0 Release，修复版 Ni/Zn 各 88 个轨道通过有限性、归一化及首次
h/i 节点检查，SETPOT/TAIL 警告均消失。Z=28、Ne=4 的旧对照轨道能量保持一致。
当时最终 CTest 为 10/10；曾出现一次 `rmcdhf_mpi` 退出阶段段错误，新旧初猜文件逐字节
一致，旧版对照和完整重跑均通过。该偶发退出问题未由本次修复解决。

初次迁移提交 `1edcfd6` 使用 GNU Fortran **15.3.0**、本机 Open MPI **5.0.11**，在独立 Release/Check
目录验证。Check 开启数组运行时检查以及 invalid/zero/overflow 浮点异常 trap。
结果如下（不包括 CSFG 扩展）：

| 验证 | 结果 |
| --- | --- |
| GNU 15 Release 全量编译/链接 | 通过 |
| GNU 15 Check 的 rwfnestimate 编译/链接 | 通过 |
| Release 专项 CTest | 3/3 通过 |
| Check 专项 CTest | 3/3 通过，含运行时检查和浮点异常 trap |
| header 输入脚本 | 9 个测试方法，含三种初猜方式、有效/无效输入及多个子用例，均通过 |
| Release 全量 CTest | 81 项中 77 项通过、4 项 Linux 专用测试跳过；其中一项路径测试在规范化 TMPDIR 后复核通过 |
| 原生 rCSFs V2 接口对照 | 32 个 CSFs、Z=28/Ne=4；完整文本/header 初猜逐字节一致，header 与描述符文件哈希均不变 |

全量首次运行中 `f2018_matrix_warnings_driver` 因 macOS 临时目录的 `/var` 与
`/private/var` 别名比较失败；该测试及警告采集逻辑未被本次修改。
使用 `TMPDIR=/private/tmp` 单独重跑此项通过，未为此改动通用测试代码。
`p0_trace_mpi_1/2`、`strict_mem_restart_mpi_1/2` 明确要求 Linux `/proc/meminfo`，
本机依设计跳过，不能算作已验证。

原生 rCSFs 对照的 header SHA-256 为
`e27c5c7fa694dbd5c690ebe5a2aeda239699cc5612f250863a00da2429949a98`，
与描述符 Parquet footer 的 `source_header_sha256` 一致；suite 读取前后保持不变。
该对照使用上游 example9 的公开小型输入，未把用户的大型计算数据复制进仓库。

可复现的专项入口：

```fish
ctest --test-dir build-rwfn-header-release -R rwfnestimate --output-on-failure
ctest --test-dir build-rwfn-header-check -R rwfnestimate --output-on-failure
```

核数据与轨道生成器在 [test/rwfnestimate](../../test/rwfnestimate)，避免依赖临时目录或
外部计算数据。测试打印每次运行的输入/日志/波函数目录；CTest 总日志在构建目录的
`Testing/Temporary/LastTest.log`。高 l 回归脚本的 `--workdir` 参数可指定长期保存的位置。
本轮构建、全量/专项 CTest 和原生接口对照记录另保存在本地
`build-rwfn-header-release/validation/`，这些运行产物不纳入 Git。

### 2026-11-12 远程同步复核

合并本地 `1edcfd6` 与远程 `e82abc0`（含 `212bcc6` 的目录重命名和
`db39505` 的 strict/diag 退役）。初猜程序及 TOML 读取模块迁至
`src/appl/rwfnestimate/`，实现内容与本地提交逐字节一致；文档索引、源码链接和
测试清单保留双方更新。前述 81 项记录对应初次迁移版本，已退役的测试不属于当前清单。

使用同一 GNU Fortran 15.3.0/Open MPI 5.0.11 工具链，在全新
`build-sync-20261112/` 中完成普通版 Release 全量编译，启用串行、MPI 和两个
OpenMP 目标，CSFG 关闭。F2018 目标/测试清单检查通过；设置
`TMPDIR=/private/tmp` 后，全量 CTest **66/66 通过，无失败或跳过**，包括 Ni/Zn、
TOML 输入和远程新增的 MPI 描述符回归。构建及测试日志保存在该构建目录的
`validation/` 和 `Testing/Temporary/LastTest.log`，不纳入 Git。

上述结论针对初猜生成和输入一致性。下一步科学验证应使用用户的真实 Ni/Zn CSF 模型，
比较后续 MCDHF 优化收敛、网格敏感性和最终能量；本记录不宣称这些物理检查已完成。
