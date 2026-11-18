# rlevels 能级表

`rlevels` 从 RMCDHF 的 `name.m` 或 RCI 的 `name.cm` 读取分块 `G92MIX` 文件。
Suite 默认构建这个独立的 Fortran 程序，安装位置为 `bin/rlevels`。
在已配置的构建目录中，可以单独构建、测试和安装该程序。

```fish
cmake --build build --target rlevels
ctest --test-dir build -R '^rlevels_file_output$' --output-on-failure
cmake --install build --component rlevels
```

```fish
./build/bin/rlevels name.m > name.level
./build/bin/rlevels first.cm second.cm > combined.level
```

省略命令行参数时，可以从标准输入逐行提供文件名，输入空行或到达 EOF 后输出能级表。
能级表写到标准输出，诊断写到标准错误。混合文件无法读取、结构不符或输入总计没有能级时，程序返回非零退出状态。

## 与上游的兼容性

实现基于上游 GRASP `488e875d` 的 `src/tool/rlevels.f90`，原程序由 Xinghong He 编写，G. Gaigalas 增加 LSJ 标记。
程序读取每个对称性块的根编号及本征值，将块能量偏移加回本征值，并消耗混合系数记录。
RCI 文件可以包含零根块头，RMCDHF 文件可以省略零根块。所有输入的总能量按升序合并，简并能级按输入顺序排列。

`Pos` 是混合文件中的块内根编号。`Energy Total` 以 Hartree 输出。`Levels` 是相对于合并表中最低能级的激发能，
`Splitting` 是相对于前一行的能差，两者沿用上游的换算因子 `2 × 109737.31568508 cm⁻¹/Hartree`。
最低能级的两列都为零，没有构型标签时也可直接读取这两列。

可选的 `name.lsj.lbl` 必须与 `.m` 或 `.cm` 同目录、同文件名主体。
构型列取 `jj2lsj` 或 `jj2lsj_omp` 所列的第一条构型，仅在对应混合文件内匹配块内序号、`J`、宇称及带符号总能量，容差为 `10⁻⁸` Hartree。
某个能级缺少匹配标签时，构型列留空。标签位置字段为 `***` 时，该行没有可用于匹配的序号。

CTest 的 `rlevels_file_output` 使用仓库中的 `data/test.m` 与 `data/test.lsj.lbl` 校验能级文本，
并生成临时 `.m` 与 `.cm` 文件检查多文件排序、空块、简并标签、文件名和错误诊断。
