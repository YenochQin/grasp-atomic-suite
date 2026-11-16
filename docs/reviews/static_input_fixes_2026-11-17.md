# RIS4 缓存与 JJ2LSJ 唯一标签输入保护

日期：2026-11-17。源码基准：`3e02e0f407d085504d9b2fb8445d50e8d013b37c`。
本记录对应[静态诊断计划](../plans/static_diagnostics_fixes_2026-11-16.md)的 S1/S2 首批代码与过程级回归。
S3 数值常量、S4 接口与防御性工作仍待独立实施；CSF reader 归描述符迁移任务。

后续进展：S3/S4 已按独立批次实施，见[精度与接口记录](static_numeric_interfaces_2026-11-17.md)。
本报告保留 S1/S2 交付时的范围与验证结果。

## 行为与根因

RIS4 原缓存探测忽略 BACKSPACE 状态，尾记录 READ 没有 IOSTAT，空文件会触发运行时 EOF。
现在缺文件、空文件、截断/短记录及非结束尾记录返回不可用，并关闭连接，允许重新计算。
打开或原始读取失败、有效尾记录后的 REWIND 失败、CLOSE 失败通过 status 返回。
额外的一字节 stream 读取区分基础可读性错误与 sequential record 结构不完整；例如目录可能 OPEN 成功而 READ 失败。
尾部结构探测失败仍按不完整缓存处理，不承诺检测所有底层文件系统错误或验证整个缓存内容。
有效缓存保持原约定：连接在 unit 50/51，已回卷，可供调用者读取。
串行调用者使用现有 `fail_input_io`；MPI root 使用 `fail_input_io_mpi`，失败时 MPI_Abort，
正常时仍在广播可用性之后使用结果。root 关闭有效缓存的状态也被检查。

`uniquelsj` 原来未检查头/能级记录，只把混合系数读取的 -1 当作 EOF，正转换错误后继续使用字段。
现在使用命名 `iostat_end`，每次成功读整行后先识别下一能级或空白组分隔，再按相应格式转换。
这避免将合法的下一能级行误当作坏分量。空文件、无能级、空分量、短记录、转换错误、
非有限值、负权重、无正的候选权重和标签耗尽均带 `uniquelsj:` 上下文非零退出。
仅在完整分量之后的 EOF 正常结束；空白行开始下一组，缺失该组能级仍是错误。
零值分量不再冒充空白分隔记录；冲突搜索限制在实际读入的分量内，防止访问未初始化标签。
唯一标签选择/重排算法、输出格式及 10000×10000 生产容量保留，正常返回释放本次工作区。

## 验证方法与结果

- 基线生产 `angdata` 在空缓存上以运行时 EOF 失败，缺失和合法缓存的控制通过。
- 基线完整 `uniquelsj` 在空输入上产生缺少过程上下文的运行时失败。
  私有过程逐字提取到测试 host，保留整个过程和生产容量，并链接实际 profile/input_io 模块；
  不将单条 READ 探针称为完整过程测试，也没有公开私有入口来方便测试。
- 三个正常唯一标签夹具（单能级、标签冲突重分配、多组）两份输出与基线逐字节一致。
  六份基线输出固定保存在 `test/jj2lsj/unique_expected/`，CTest 持续执行字节比较。
  另有 11 个损坏/不可分配标签的输入场景检查非零状态及错误上下文；全程序正常路径由既有
  `jj2lsj_parallel_regression` 覆盖。
- RIS4 每种串行/MPI 驱动覆盖 16 个 one-body/two-body 场景：缺失、结束记录、
  带正常系数的缓存、非结束记录、空、两字节截断、短记录、目录读取错误。
  可恢复场景连续调用两次，检查不可用缓存已关闭、有效缓存已回卷、系数 `0.125` 保持不变。
  MPI 使用两个 rank，仅 root 读文件，随后广播；root 读取失败要求整组在 20 秒内非零退出且有上下文。
  该驱动调用实际 `angdata` 和实际 MPI 错误处理模块，并非完整 `ris_cal_mpi` 原子计算。
- GNU Fortran 15.3.0、Open MPI 5.0.11、Fortitude 0.10.1。本机规范指定的
  `/usr/bin/gfortran-15` 不存在，使用已安装 Homebrew GNU 15，并显式设置 FC/OMPI_FC，
  未修改仓库工具链规范。两套新目录均全量构建，未安装覆盖生产二进制。
- Check 首轮 79/81；两项既有测试因 macOS 临时目录 `/var` 与 `/private/var` 路径别名失败。
  指定规范化 TMPDIR 后，两项和扩展后的四项新增测试共 6/6 通过。
- Release 全量 81/81 通过。最后清理唯一标签错误提示后，两套构建的 JJ2LSJ 五项均复测通过。
- 四份修改的生产实现定向 C181 检查通过。C183 仍报告缺少 IOMSG；当前已有阶段、unit、
  状态码及可获取的文件名诊断，本轮没有以清空该建议为目标。没有重新运行全库 fortls 或放宽告警预算。

未新增完整 RIS4 原子物理结果对照，未注入难以稳定复现的文件系统 REWIND/CLOSE 故障。
缓存内部所有记录的语义正确性、生产规模性能和物理收敛不是本轮已验证结论。

## 复查入口

新增测试位于 `test/ris4/{cache_io.F90,check_cache_io.py}` 和
`test/jj2lsj/{unique_io.f90,check_unique_io.py}`；CMake 与 `test/f2018/required.json` 已登记。

```fish
ctest --test-dir build -R 'ris4_cache|jj2lsj_unique|jj2lsj_parallel' --output-on-failure
```

本机忽略目录 `build-static-input-20261117/` 保存基线源码、驱动、正常输出、失败日志，
以及 `check-build.log`、`release-build.log`、`check-ctest.log`、`check-retest.log`、
`release-ctest.log` 和 `c181.log`。这些不是仓库可移植依赖。
