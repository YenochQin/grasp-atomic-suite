# grasp-atomic-suite

`grasp-atomic-suite` 是一个基于 GRASP2018/RHFS 代码体系扩展的研究型工作仓库，面向相对论原子结构与原子性质计算的整理、分析、验证与开发。

这个仓库不是单纯的上游镜像。除了保留现有 Fortran 源码树外，还覆盖了 `rhfs90`、`ris4`、`gj90` 以及 MPI 并行相关组件，并保留围绕 `Landé g_J` 计算链路的分析文档与后续拆分设计。

## 当前内容

- 保留了 GRASP 风格的源代码结构，便于直接复用现有数值核与构建流程
- 包含 `rhfs90`、`ris4`、`gj90` 等程序目录及其相关依赖
- 包含 `mpi90` 等并行支持模块，用于 RHFS 和 RIS 的并行计算
- 分析了 `rhfs90` 中 `g_J` 的实现链路
- 记录了 `g_J` 的理论背景、MCDHF/RCI 下的计算思路和使用要点
- 给出了把 `g_J` 逻辑从 `rhfs90` 中拆分为独立程序 `gj90` 的设计方案

## 仓库定位

如果你只是想使用上游通用版 GRASP，这个仓库并不是最简入口。  
如果你关注以下问题，这个仓库更合适：

- 如何在一个工作仓库中组织 `rhfs90`、`ris4` 与相关并行程序
- 如何在 GRASP/RHFS 代码树基础上扩展新的原子性质计算或独立程序
- RHFS 与 RIS 的串行/并行构建入口、源码组织和依赖关系
- `rhfs90` 里 `g_J` 是如何计算出来的
- `.c`、`.w`、`.m/.cm` 文件在 `g_J` 计算中的具体作用
- 如何从 RHFS 的现有实现中抽出一个只计算 `g_J` 的程序
- 如何在后续开发中验证 `gj90` 与 `rhfs90` 的数值一致性

## 目录概览

```text
.
├── README.md
├── lande_g.md
├── RHFS_gJ_report.md
├── gj90_blueprint.md
├── configure.sh
├── grasptest/
├── src/
│   ├── appl/
│   │   ├── gj90/
│   │   ├── rdensity/
│   │   ├── rhfs90/
│   │   ├── rhfszeeman95/
│   │   └── ris4/
│   └── lib/
│       ├── mpi90/
│       └── ...
├── bin/
└── lib/
```

## 重点文件

- `lande_g.md`
  从理论角度说明 `Landé g_J` 的定义、MCDHF/RCI 中的实现思路、与 LS 近似公式的关系，以及在 GRASP 中的实际使用注意事项。

- `RHFS_gJ_report.md`
  追踪 `src/appl/rhfs90` 中 `g_J` 的具体实现路径，包含入口程序、输入文件读取、径向积分、角向矩阵元与 ASF 投影过程。

- `gj90_blueprint.md`
  给出独立程序 `gj90` 的设计目标、最小依赖链、建议复用模块以及分阶段开发计划。

- `src/appl/rhfs90/`
  当前 `g_J` 计算相关源码的核心目录，重点文件包括：
  - `hfs92.f90`
  - `hfsgg.f90`
  - `gethfd.f90`
  - `getmixblock.f90`
  - `matelt.f90`
  - `rinthf.f90`

## 构建方法

推荐使用 CMake 的 out-of-source 构建：

```sh
./configure.sh
cd build
make -j4 install
```

调试构建：

```sh
./configure.sh --debug
cd build-debug
make -j4 install
```

清理构建产物：

```sh
./scripts/clean-build-artifacts.sh
```

也保留了传统 `Makefile` 工作流：

```sh
make
make src/lib/libmod
make src/appl/rhfs90
```

说明：

- CMake 构建产物会先进入 `build/` 或 `build-debug/`，执行 `make install` 后安装到仓库根目录下的 `bin/` 和 `lib/`
- 传统 `Makefile` 会直接在源码树附近生成中间文件，并把程序安装到 `bin/`、`lib/`
- 本仓库已有 `.gitignore`，默认忽略常见编译产物和构建目录

## 编译配置

不要直接修改顶层 `Makefile`。建议复制模板生成 `Make.user`：

```sh
cp Make.user.gfortran Make.user
```

或者：

```sh
cp Make.user.ifort Make.user
```

然后按本机环境调整：

- `FC`
- `FC_MPI`
- `FC_FLAGS`
- `FC_LD`

## 与 `g_J` 相关的输入文件

围绕 `rhfs90` 的 `g_J` 计算，核心输入通常包括：

- `name.c`
  CSF 列表与耦合信息

- `name.w`
  径向轨道波函数

- `name.m` 或 `name.cm`
  ASF 混合系数与本征信息

- `isodata`
  同位素和径向网格相关输入

这些文件的读取和作用分工在 `RHFS_gJ_report.md` 中有更详细的逐步说明。

## 测试与示例

示例和回归式脚本位于 `grasptest/`。常见用法：

```sh
grasptest/example1/script/script_ex1
grasptest/case1/script/sh_case1
```

如果本地构建目录中包含 `test/`，也可以运行：

```sh
cd build
ctest
```

## 当前状态

当前仓库已经完成：

- 保留并整理了面向 `rhfs90`、`ris4` 和并行构建的源码树
- `rhfs90` 中 `g_J` 计算链路的源码梳理
- `g_J` 理论与程序实现之间的对应关系说明
- 独立程序 `gj90` 的初步设计

当前仓库尚未承诺：

- `gj90` 已经实现并可直接生产使用
- 所有文档都与上游 GRASP 各版本完全同步
- 对所有原子体系都给出完整验证数据

## 参考背景

本仓库基于 GRASP2018 风格代码树展开工作。GRASP2018 相关论文：

> C. Froese Fischer, G. Gaigalas, P. Jönsson, J. Bieroń,  
> "GRASP2018 — a Fortran 95 version of the General Relativistic Atomic Structure Package",  
> Computer Physics Communications, 237, 184-187 (2018),  
> https://doi.org/10.1016/j.cpc.2018.10.032

## 后续方向

- 继续整理 `rhfs90`、`ris4` 及其 MPI 版本的构建和计算链路
- 从 `rhfs90` 中分离只计算 `g_J` 的数值核心
- 建立 `gj90` 的最小可运行版本
- 用现有 `rhfs90` 输出逐态核对 `g_J`、`delta g_J` 和 `total g_J`
- 补充更多围绕原子结构与原子性质计算的说明、测试与示例

## 许可证

仓库代码遵循 [MIT License](LICENSE)。
