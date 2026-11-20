# spin_angular 变量与常量命名迁移

日期 2026-11-20。迁移基准为 suite 的 `2d4c721` 与当前工作树中的 CSF 查询合并及共享 kind
导入清理。盘点原 `librang90` 中参与构建的 **89 份 Fortran 源码**。
文件与模块现名见[自旋–角向库命名迁移](spin_angular_naming.md)。
每个过程按 Fortran 不区分大小写的规则登记原名，按作用域区分同名实体。

- 本地声明共 **1,825 项**，重命名 **1,819 项**，沿用 **6 项**。
  范围含 660 个形参（含 1 个过程形参）、1,073 个局部量、
  82 个 DATA 初始化对象和 10 个函数结果。此目录没有本地 `PARAMETER` 声明。
- `USE ONLY` 按作用域和局部名称去重后，共登记 **784 处**数据和 kind 导入，其中 **510 处**使用局部别名。
  表列出各别名的定义模块和使用过程。共享数组继续由相应模块持有，`real64`、`max_orbitals` 和
  `zero`、`half` 等已有清晰名称的常量式用量沿用现名。DATA 对象保留原初始化、数组范围和存储属性。
- 例程名、参数顺序和传递方式沿用原接口。对外部 Fortran 源码，位置实参可沿用，
  关键字实参需要采用下表中的新名称。数字字面量、文件标签和诊断字符串沿用原值。
- `CMakeLists.txt` 声明 `spin_angular` 静态库，没有本地 CMake 变量定义。`ReadMe` 是历史说明，
  其中的标识符用于追溯旧代码。

## 命名依据

- **CSF、轨道、子壳层位置。** `bra_csf_index` 和 `ket_csf_index` 是左右 CSF 编号。
  `*_orbital_index` 及 `orbital_a` 等索引完整轨道表，`*_shell_position` 及 `bra_shell_1`
  等索引活跃子壳层列表。`peel_orbital_indices` 将子壳层位置转换为完整轨道号，
  临时插入闭壳层时它也承载插入后的次序。
  见 [RKCO_GG](../../src/lib/spin_angular/two_electron_angular_coefficients.f90) 和 [ONEPARTICLEJJ](../../src/lib/spin_angular/one_electron_tensor_coefficients.f90)。
- **角动量与准自旋编码。** `two_orbital_j` 表示单电子的 `2j = 2|κ| − 1`，
  `two_term_j` 表示子壳层项的 `2J`，`two_q` 表示准自旋的 `2Q`。
  `two_j_plus_one` 和各 `*_j_weight` 保留 `2j + 1` 或 `2J + 1` 的权重编码。
  `m_c::jjc1/jjc2` 中的中间耦合采用 `2J + 1`。见 [DIAGA2](../../src/lib/spin_angular/later_active_shell_recoupling.f90)。
  Wigner 6j 的 `two_j_1` 到 `two_j_6` 按参数位置编号。Wigner 9j 的 `two_j_23` 表示第 2 行、第 3 列。
- **子壳层张量数据。** `*_shell_codes` 的 7 槽依次为项索引、主量子数、`2j`、占据数、
  轨道角动量 `l`、子壳层 `2J`、准自旋 `2Q`。`*_shell_values` 的 3 槽依次为
  `Q`、`J`、准自旋投影 `M_Q`。`q_projection_*` 是算符的准自旋投影。
  `mt*` 将 `2Q`、两倍 seniority 与 `2J` 打包为整数。
  见 [PERKO1](../../src/lib/spin_angular/jj_subshell_state_packing.f90)、[NMTEJJ](../../src/lib/spin_angular/jj_subshell_term_lookup.f90)
  和 [RUMTJJ](../../src/lib/spin_angular/jj_subshell_term_decoding.f90)。
- **算符、系数与积分编号。** `A` 表示单个产生或湮灭张量，`W` 表示两个 `A` 耦合而成的张量。
  `cg_coefficient` 为 Clebsch–Gordan 系数，`reduced_cfp` 为完全约化的分数亲代系数。
  `interaction_type=1` 走 Coulomb 路径，`interaction_type=2` 走 Breit 路径。
  Breit 的 `first_direct_nu`、`first_exchange_nu` 分别为两组径向阶数起点，
  `direct_nu_count`、`exchange_nu_count` 是步长为 2 的候选个数。
  见 [WJ1](../../src/lib/spin_angular/jj_w_matrix_element.f90)、[RMEAJJ](../../src/lib/spin_angular/jj_reduced_cfp.f90)
  和 [SNRC](../../src/lib/spin_angular/breit_radial_rank_ranges.f90)。
- **选择规则与求值。** 重耦合例程的整型 `evaluate_coefficient` 为 0 时检查允许的耦合，
  非零时求重耦合因子。若原变量在多个阶段复用，下表说明其不同角色。仅声明的历史局部量照常登记。

## 验证

- 将新标识符和局部导入别名逆映射，并把显式 RESULT 归一为原隐式结果绑定。
  89 份源码逐语句的 Fortran 标记序列与迁移前相同。检查了声明覆盖、
  同名冲突、续行、字符串、DATA 初始化及 OpenMP THREADPRIVATE 子句。扫描仓库 `src/`
  和 `test/` 的 525 处相关过程引用，均使用位置实参。
- GNU Fortran 15.3.0 构建通过。启用 MPI 和 jj2lsj OpenMP 的 Release 配置通过 110 项 CTest。
  关闭 MPI 和 OpenMP 的 Check 配置通过 70 项 CTest，Rust reader 测试在 Release 配置中执行。
- 使用同配置的迁移前后 Release 二进制对照相同输入夹具。
  串行、2 rank 和 4 rank 的 `rangular` 分别有 6、12、24 个 MCP 文件逐字节一致。
  `rci_mpi` 的 2 rank Ni I 七 CSF 夹具分别在 Coulomb 和 Breit 路径运行，
  两条路径各自的 `case.cm` 均与迁移前逐字节一致。

## 共享状态与导入常量

同一模块实体在不同过程中的局部名称一致。`ribojj*_c` 中的 `same_parity` 和
`opposite_parity` 是相对当前子壳层占据配对使用的项索引区间。

| 定义模块 | 原局部名 | 迁移后 | 使用过程 |
| --- | --- | --- | --- |
| `cons_c` | `eleven` | `eleven` | `oneparticlejj`, `onescalar`, `rwjj` |
| `cons_c` | `eps` | `angular_tolerance` | `a1jj`, `awp1`, `awp1jjg`, `c0t5s`, `c1e0sm`, `c1e1sm`, `cle0sm`, `coulom`, `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1112`, `gg1122`, `gg12`, `gg1222`, `gg1233`, `gg1234`, `oneparticlejj`, `oneparticlejj1`, `oneparticlejj2`, `onescalar`, `onescalar1`, `onescalar2`, `rmeajj11`, `rwjj`, `sixj4`, `suwjj`, `w1jjg`, `wap1`, `wap1jjg`, `wj1`, `ww1` |
| `cons_c` | `four` | `four` | `oneparticlejj`, `onescalar`, `rwjj`, `sixj2`, `sixj3`, `sixj35` |
| `cons_c` | `half` | `half` | `c0t5s`, `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1112`, `gg1222`, `itjj`, `itjj2`, `itjj3`, `oneparticlejj`, `oneparticlejj1`, `oneparticlejj2`, `onescalar`, `onescalar1`, `onescalar2`, `perko1`, `rmeajj11`, `rwjj`, `sixj1`, `sixj2`, `sixj4`, `w1jjg` |
| `cons_c` | `one` | `one` | `c0t5s`, `c1e0sm`, `c1e1sm`, `cle0sm`, `diaga3`, `gracah1`, `nine0`, `oneparticlejj`, `onescalar`, `onescalar1`, `reco2`, `reco3`, `reco4`, `recop1`, `recop2`, `rwjj`, `sixj`, `sixj1`, `sixj2`, `sixj3`, `sixj35`, `sixj4`, `sixj5`, `wj1` |
| `cons_c` | `seven` | `seven` | `oneparticlejj`, `onescalar`, `rwjj`, `sixj3` |
| `cons_c` | `ten` | `ten` | `oneparticlejj`, `onescalar`, `rwjj` |
| `cons_c` | `tenth` | `tenth` | `awp1`, `awp1jjg`, `c0t5s`, `c1e0sm`, `c1e1sm`, `cle0sm`, `gg1112`, `gg1122`, `gg12`, `gg1222`, `gg1233`, `gg1234`, `izas1`, `oneparticlejj`, `oneparticlejj2`, `onescalar`, `rwjj`, `w1jjg`, `wap1`, `wap1jjg`, `wj1`, `ww1` |
| `cons_c` | `three` | `three` | `c1e0sm`, `c1e1sm`, `oneparticlejj`, `onescalar`, `rwjj`, `sixj1`, `sixj2`, `sixj3`, `sixj35`, `sixj4` |
| `cons_c` | `two` | `two` | `awp1`, `awp1jjg`, `c0t5s`, `c1e0sm`, `c1e1sm`, `cle0sm`, `gracah1`, `izas1`, `oneparticlejj`, `onescalar`, `rwjj`, `sixj1`, `sixj2`, `sixj3`, `sixj35`, `sixj4`, `sixj5`, `w1jjg`, `wap1`, `wap1jjg`, `wj1`, `ww1` |
| `cons_c` | `zero` | `zero` | `a1jj`, `awp1`, `awp1jjg`, `c0t5s`, `c1e0sm`, `c1e1sm`, `cle0sm`, `coulom`, `diaga1`, `diaga2`, `diaga3`, `diaga4`, `diaga5`, `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1112`, `gg1122`, `gg12`, `gg1222`, `gg1233`, `gg1234`, `gracah1`, `nine`, `nine0`, `oneparticlejj`, `oneparticlejj1`, `oneparticlejj2`, `onescalar`, `onescalar1`, `onescalar2`, `rmeajj`, `rmeajj11`, `rmeajj9`, `rmew1jj`, `rmew3jj`, `rmew5jj`, `rmew7bjj`, `rmew7jj`, `rwjj`, `sixj`, `sixj1`, `sixj2`, `sixj3`, `sixj35`, `sixj4`, `sixj5`, `suwjj`, `w1jjg`, `wap1`, `wap1jjg`, `wj1`, `ww1` |
| `dumx_c` | `jc1s` | `saved_bra_j_weights` | `oneparticlejj`, `onescalar`, `rkco_gg` |
| `dumx_c` | `jc2s` | `saved_ket_j_weights` | `oneparticlejj`, `onescalar`, `rkco_gg` |
| `dumx_c` | `jlis` | `saved_peel_orbitals` | `oneparticlejj`, `onescalar`, `rkco_gg` |
| `facts_c` | `gam` | `log_factorials` | `dracah`, `gracah1` |
| `facts_c` | `mfact` | `factorial_table_capacity` | `dracah`, `gracah1` |
| `iso_fortran_env` | `real64` | `real64` | `a1jj`, `awp1`, `awp1jjg`, `c0t5s`, `c1e0sm`, `c1e1sm`, `cle0sm`, `coulom`, `diaga1`, `diaga2`, `diaga3`, `diaga4`, `diaga5`, `dracah`, `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1112`, `gg1122`, `gg12`, `gg1222`, `gg1233`, `gg1234`, `gracah1`, `itjj`, `itjj2`, `itjj3`, `izas1`, `nine`, `nine0`, `oneparticlejj`, `oneparticlejj1`, `oneparticlejj2`, `onescalar`, `onescalar1`, `onescalar2`, `perko1`, `rec3`, `reco2`, `reco3`, `reco4`, `recop1`, `recop2`, `rmeajj`, `rmeajj11`, `rmeajj9`, `rmew1jj`, `rmew3jj`, `rmew5jj`, `rmew7bjj`, `rmew7jj`, `rwjj`, `sixj`, `sixj1`, `sixj2`, `sixj3`, `sixj35`, `sixj4`, `sixj5`, `suwjj`, `w1jjg`, `wap1`, `wap1jjg`, `wj1`, `ww1` |
| `m_c` | `jjc1` | `bra_intermediate_couplings` | `diaga1`, `diaga2`, `diaga3`, `diaga4`, `diaga5`, `oneparticlejj`, `onescalar`, `onescalar2`, `reco`, `recoonescalar`, `recop00`, `rkco_gg` |
| `m_c` | `jjc2` | `ket_intermediate_couplings` | `diaga1`, `diaga2`, `diaga3`, `diaga4`, `diaga5`, `oneparticlejj`, `onescalar`, `onescalar2`, `reco`, `recoonescalar`, `recop00`, `rkco_gg` |
| `m_c` | `jjq1` | `bra_subshell_quantum_codes` | `diaga1`, `diaga2`, `diaga3`, `diaga4`, `diaga5`, `oneparticlejj`, `onescalar`, `onescalar2`, `perko1`, `reco`, `reco2`, `reco3`, `reco4`, `recoonescalar`, `recop00`, `recop1`, `recop2`, `rkco_gg` |
| `m_c` | `jjq2` | `ket_subshell_quantum_codes` | `diaga1`, `diaga2`, `diaga4`, `diaga5`, `oneparticlejj`, `onescalar`, `onescalar2`, `perko1`, `reco`, `reco2`, `reco3`, `reco4`, `recoonescalar`, `recop00`, `recop2`, `rkco_gg` |
| `m_c` | `jlist` | `peel_orbital_indices` | `diaga1`, `diaga2`, `diaga3`, `diaga4`, `diaga5`, `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `oneparticlejj`, `oneparticlejj2`, `onescalar`, `onescalar2`, `perko1`, `reco`, `reco2`, `reco3`, `reco4`, `recoonescalar`, `recop00`, `recop1`, `recop2`, `rkco_gg` |
| `m_c` | `klist` | `core_orbital_indices` | `oneparticlejj`, `onescalar`, `onescalar2`, `rkco_gg` |
| `m_c` | `max_orbitals` | `max_orbitals` | `onescalar2` |
| `m_c` | `ncore` | `core_orbital_count` | `oneparticlejj`, `onescalar`, `onescalar2`, `rkco_gg` |
| `m_c` | `npeel` | `peel_orbital_count` | `el1`, `el2`, `el3`, `el31`, `el32`, `el33`, `el4`, `el41`, `el5`, `el51`, `el52`, `el53`, `oneparticlejj`, `onescalar`, `onescalar2`, `reco`, `recoonescalar`, `recop00`, `recop1`, `recop2`, `rkco_gg` |
| `m_c` | `nq1` | `bra_occupations` | `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `oneparticlejj`, `oneparticlejj2`, `onescalar`, `onescalar2`, `perko1`, `rkco_gg` |
| `m_c` | `nq2` | `ket_occupations` | `oneparticlejj`, `onescalar`, `onescalar2`, `perko1`, `rkco_gg` |
| `m_c`, `orb_c`, `dumx_c` | `max_orbitals` | `max_orbitals` | `rkco_gg` |
| `mtjj2_c` | `mt11` | `packed_term_codes_high_j` | `nmtejj`, `rumtjj` |
| `mtjj2_c` | `mt9` | `packed_term_codes_two_j_9` | `nmtejj`, `rumtjj` |
| `mtjj_c` | `mt` | `packed_term_codes_low_j` | `nmtejj`, `rumtjj` |
| `orb_c` | `orbital_count` | `orbital_count` | `oneparticlejj`, `onescalar`, `rkco_gg` |
| `orb_c` | `orbital_kappa` | `orbital_kappa` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `oneparticlejj`, `oneparticlejj2`, `onescalar`, `perko1`, `rkco_gg` |
| `orb_c` | `orbital_principal_n` | `orbital_principal_n` | `perko1` |
| `parameter_def`, `m_c` | `max_orbitals` | `max_orbitals` | `oneparticlejj`, `onescalar` |
| `ribojj11_c` | `imgnjj11` | `last_opposite_parity_terms_high_j` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj11_c` | `imgtjj11` | `last_same_parity_terms_high_j` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj11_c` | `impnjj11` | `first_opposite_parity_terms_high_j` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj11_c` | `imptjj11` | `first_same_parity_terms_high_j` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj9_c` | `imgnjj9` | `last_opposite_parity_terms_two_j_9` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj9_c` | `imgtjj9` | `last_same_parity_terms_two_j_9` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj9_c` | `impnjj9` | `first_opposite_parity_terms_two_j_9` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj9_c` | `imptjj9` | `first_same_parity_terms_two_j_9` | `itjj`, `itjj2`, `itjj3`, `w1jjg` |
| `ribojj_c` | `imgnjj` | `last_opposite_parity_terms_low_j` | `itjj`, `itjj2`, `itjj3`, `rmeajj`, `rmew1jj`, `rmew3jj`, `rmew5jj`, `rmew7bjj`, `rmew7jj`, `rwjj`, `suwjj`, `w1jjg` |
| `ribojj_c` | `imgtjj` | `last_same_parity_terms_low_j` | `itjj`, `itjj2`, `itjj3`, `rmeajj`, `rmew1jj`, `rmew3jj`, `rmew5jj`, `rmew7bjj`, `rmew7jj`, `rwjj`, `suwjj`, `w1jjg` |
| `ribojj_c` | `impnjj` | `first_opposite_parity_terms_low_j` | `itjj`, `itjj2`, `itjj3`, `rmeajj`, `rmew1jj`, `rmew3jj`, `rmew5jj`, `rmew7bjj`, `rmew7jj`, `rwjj`, `suwjj`, `w1jjg` |
| `ribojj_c` | `imptjj` | `first_same_parity_terms_low_j` | `itjj`, `itjj2`, `itjj3`, `rmeajj`, `rmew1jj`, `rmew3jj`, `rmew5jj`, `rmew7bjj`, `rmew7jj`, `rwjj`, `suwjj`, `w1jjg` |
| `trk_c` | `bd1` | `ket_shell_values_1` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bd2` | `ket_shell_values_2` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bd3` | `ket_shell_values_3` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bd4` | `ket_shell_values_4` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bk1` | `bra_shell_values_1` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bk2` | `bra_shell_values_2` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bk3` | `bra_shell_values_3` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `bk4` | `bra_shell_values_4` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `id1` | `ket_shell_codes_1` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `id2` | `ket_shell_codes_2` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `id3` | `ket_shell_codes_3` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `id4` | `ket_shell_codes_4` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `ik1` | `bra_shell_codes_1` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `ik2` | `bra_shell_codes_2` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `ik3` | `bra_shell_codes_3` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |
| `trk_c` | `ik4` | `bra_shell_codes_4` | `el1`, `el2`, `el31`, `el32`, `el33`, `el41`, `el51`, `el52`, `el53`, `gg1122`, `oneparticlejj1`, `oneparticlejj2`, `onescalar1`, `onescalar2`, `perko2` |

## 逐过程对照

### [a1jj](../../src/lib/spin_angular/jj_a_matrix_element_few_electrons.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `qm1` | `operator_q_projection_1` | 形参 |  |
| `a` | `a_matrix_element` | 形参 |  |
| `ifaz` | `twice_phase_exponent` | 局部量 |  |
| `isuma` | `integer_norm_squared` | 局部量 |  |
| `ab` | `norm_squared` | 局部量 |  |

### [awp1](../../src/lib/spin_angular/jj_aw_matrix_element.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `bk2` | `total_angular_rank` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `aw` | `aw_matrix_element` | 形参 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `kk2` | `two_total_rank` | 局部量 |  |
| `ie` | `twice_phase_exponent` | 局部量 |  |
| `iq` | `intermediate_occupation_change` | 局部量 |  |
| `iq2` | `occupation_change_2` | 局部量 |  |
| `iq3` | `occupation_change_3` | 局部量 |  |
| `iqm` | `two_abs_intermediate_projection` | 局部量 |  |
| `it` | `intermediate_term_index` | 局部量 |  |
| `itp` | `first_intermediate_term` | 局部量 |  |
| `itg` | `last_intermediate_term` | 局部量 |  |
| `ibtt` | `two_intermediate_seniority` | 局部量 |  |
| `ibt` | `intermediate_shell_codes` | 局部量 |  |
| `enqp` | `sum_over_terms` | 局部量 |  |
| `d1` | `a_matrix_element_term` | 局部量 |  |
| `s` | `reduced_cfp` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `w` | `w_matrix_element` | 局部量 |  |
| `bt` | `intermediate_shell_values` | 局部量 |  |

### [awp1jjg](../../src/lib/spin_angular/jj_aw_matrix_element_high_j.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `bk2` | `total_angular_rank` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `aw` | `aw_matrix_element` | 形参 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `kk2` | `two_total_rank` | 局部量 |  |
| `ie` | `twice_phase_exponent` | 局部量 |  |
| `iq` | `intermediate_occupation_change` | 局部量 |  |
| `iq2` | `occupation_change_2` | 局部量 |  |
| `iq3` | `occupation_change_3` | 局部量 |  |
| `iqm` | `two_abs_intermediate_projection` | 局部量 |  |
| `it` | `intermediate_term_index` | 局部量 |  |
| `itp` | `first_intermediate_term` | 局部量 |  |
| `itg` | `last_intermediate_term` | 局部量 |  |
| `ibtt` | `two_intermediate_seniority` | 局部量 |  |
| `ibt` | `intermediate_shell_codes` | 局部量 |  |
| `enqp` | `sum_over_terms` | 局部量 |  |
| `d1` | `a_matrix_element_term` | 局部量 |  |
| `si1` | `six_j_value` | 局部量 |  |
| `w` | `w_matrix_element` | 局部量 |  |
| `bt` | `intermediate_shell_values` | 局部量 |  |

### [c0t5s](../../src/lib/spin_angular/quasispin_cg_rank_half.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `q` | `ket_q` | 形参 |  |
| `qm` | `ket_q_projection` | 形参 |  |
| `sm` | `operator_q_projection` | 形参 |  |
| `c` | `bra_q` | 形参 |  |
| `cm` | `bra_q_projection` | 形参 |  |
| `a` | `cg_coefficient` | 形参 |  |
| `iiq` | `two_ket_q` | 局部量 |  |
| `iic` | `two_bra_q` | 局部量 |  |
| `ie` | `projection_sign_index` | 局部量 |  |
| `gc` | `projection_signs` | 局部量 |  |

### [c1e0sm](../../src/lib/spin_angular/quasispin_cg_rank_one_zero_projection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `q` | `ket_q` | 形参 |  |
| `qm` | `ket_q_projection` | 形参 |  |
| `c` | `bra_q` | 形参 |  |
| `cm` | `bra_q_projection` | 形参 |  |
| `a` | `cg_coefficient` | 形参 |  |
| `iiq` | `two_ket_q` | 局部量 |  |
| `iic` | `two_bra_q` | 局部量 |  |
| `is` | `zero_projection_parity_sum` | 局部量 |  |
| `ig` | `rank_one_coupling_case` | 局部量 |  |

### [c1e1sm](../../src/lib/spin_angular/quasispin_cg_rank_one_unit_projection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `q` | `ket_q` | 形参 |  |
| `qm` | `ket_q_projection` | 形参 |  |
| `sm` | `operator_q_projection` | 形参 |  |
| `c` | `bra_q` | 形参 |  |
| `cm` | `bra_q_projection` | 形参 |  |
| `a` | `cg_coefficient` | 形参 |  |
| `ie` | `projection_sign_index` | 局部量 |  |
| `iiq` | `two_ket_q` | 局部量 |  |
| `iic` | `two_bra_q` | 局部量 |  |
| `gc` | `projection_signs` | 局部量 |  |

### [cle0sm](../../src/lib/spin_angular/quasispin_cg_zero_projection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `q` | `ket_q` | 形参 |  |
| `qm` | `ket_q_projection` | 形参 |  |
| `s` | `quasispin_rank` | 形参 |  |
| `c` | `bra_q` | 形参 |  |
| `cm` | `bra_q_projection` | 形参 |  |
| `a` | `cg_coefficient` | 形参 |  |
| `iiq` | `two_ket_q` | 局部量 |  |
| `iic` | `two_bra_q` | 局部量 |  |
| `iis` | `two_quasispin_rank` | 局部量 |  |

### [coulom](../../src/lib/spin_angular/coulomb_angular_factor.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `abs_kappa_a` | 形参 |  |
| `j2` | `abs_kappa_b` | 形参 |  |
| `j3` | `abs_kappa_c` | 形参 |  |
| `j4` | `abs_kappa_d` | 形参 |  |
| `l1` | `orbital_l_a` | 形参 |  |
| `l2` | `orbital_l_b` | 形参 |  |
| `l3` | `orbital_l_c` | 形参 |  |
| `l4` | `orbital_l_d` | 形参 |  |
| `k` | `multipole_rank` | 形参 |  |
| `aa` | `coulomb_angular_factor` | 形参 |  |
| `i` | `harmonic_rank` | 局部量 | 传给 cre 的球张量阶数，按现有整数表达式由 multipole_rank 得到。 |
| `ifaz` | `twice_phase_exponent` | 局部量 |  |

### [diaga1](../../src/lib/spin_angular/first_active_shell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ij1` | `orbital_index` | 局部量 |  |
| `ia1` | `two_bra_subshell_j` | 局部量 |  |
| `ib1` | `two_ket_subshell_j` | 局部量 |  |
| `k1` | `coupling_position` | 局部量 |  |
| `j1` | `two_spectator_j` | 局部量 |  |
| `it1` | `two_bra_coupled_j` | 局部量 |  |
| `it1s` | `two_ket_coupled_j` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ll1` | `spectator_orbital` | 局部量 |  |
| `a1` | `diagram_factor` | 局部量 |  |

### [diaga2](../../src/lib/spin_angular/later_active_shell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `ia1` | `two_bra_shell_j_1` | 局部量 |  |
| `ia2` | `two_bra_shell_j_2` | 局部量 |  |
| `ib1` | `two_ket_shell_j_1` | 局部量 |  |
| `ib2` | `two_ket_shell_j_2` | 局部量 |  |
| `n1` | `coupled_position` | 局部量 |  |
| `n2` | `preceding_coupling_position` | 局部量 |  |
| `j2` | `two_bra_coupled_j` | 局部量 |  |
| `it2` | `two_bra_prefix_j` | 局部量 |  |
| `it2s` | `two_ket_prefix_j` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `a2` | `six_j_value` | 局部量 |  |

### [diaga3](../../src/lib/spin_angular/spectator_shell_chain_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `i` | `shell_position` | 局部量 |  |
| `ll1` | `orbital_index` | 局部量 |  |
| `ji` | `two_spectator_j` | 局部量 |  |
| `kk1` | `coupled_position` | 局部量 |  |
| `kk2` | `preceding_coupling_position` | 局部量 |  |
| `iti` | `two_bra_input_j` | 局部量 |  |
| `iti1` | `two_bra_coupled_j` | 局部量 |  |
| `itis` | `two_ket_input_j` | 局部量 |  |
| `iti1s` | `two_ket_coupled_j` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `aa` | `diagram_product` | 局部量 |  |
| `a3` | `diagram_factor` | 局部量 |  |

### [diaga4](../../src/lib/spin_angular/pair_tensor_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `k1` | `two_rank_1` | 形参 |  |
| `k2` | `two_rank_2` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `ia1` | `two_bra_shell_j_1` | 局部量 |  |
| `ia2` | `two_bra_shell_j_2` | 局部量 |  |
| `ib1` | `two_ket_shell_j_1` | 局部量 |  |
| `ib2` | `two_ket_shell_j_2` | 局部量 |  |
| `n1` | `coupled_position` | 局部量 |  |
| `n2` | `preceding_coupling_position` | 局部量 |  |
| `j2` | `two_bra_coupled_j` | 局部量 |  |
| `j2s` | `two_ket_coupled_j` | 局部量 |  |
| `it2` | `two_bra_prefix_j` | 局部量 |  |
| `it2s` | `two_ket_prefix_j` | 局部量 |  |
| `a2` | `nine_j_value` | 局部量 |  |

### [diaga5](../../src/lib/spin_angular/final_shell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `npeelgg` | `coupled_shell_count` | 形参 |  |
| `ja1` | `first_shell_position` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `iti1` | `two_bra_coupled_j` | 局部量 |  |
| `iti1s` | `two_ket_coupled_j` | 局部量 |  |
| `ij1` | `last_orbital` | 局部量 |  |
| `iti` | `two_bra_input_j` | 局部量 |  |
| `itis` | `two_ket_input_j` | 局部量 |  |
| `ji` | `two_spectator_j` | 局部量 |  |
| `a3` | `diagram_factor` | 局部量 |  |

### [dracah](../../src/lib/spin_angular/racah_coefficients.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `two_j_1` | 形参 |  |
| `j` | `two_j_2` | 形参 |  |
| `k` | `two_j_3` | 形参 |  |
| `l` | `two_j_4` | 形参 |  |
| `m` | `two_j_5` | 形参 |  |
| `n` | `two_j_6` | 形参 |  |
| `rac` | `racah_coefficient` | 形参 |  |
| `j1` | `triangle_sum_1` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j2` | `triangle_sum_2` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j3` | `triangle_sum_3` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j4` | `triangle_sum_4` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j5` | `quadrangle_sum_1` | 局部量 |  |
| `j6` | `quadrangle_sum_2` | 局部量 |  |
| `j7` | `quadrangle_sum_3` | 局部量 |  |
| `numin` | `sum_start` | 局部量 |  |
| `numax` | `sum_end` | 局部量 |  |
| `icount` | `terms_evaluated` | 局部量 |  |
| `kk` | `term_iteration` | 局部量 |  |
| `ki` | `descending_term_index` | 局部量 |  |

### [eile](../../src/lib/spin_angular/three_subshell_ordering.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jc` | `shell_position_c` | 形参 |  |
| `jaa` | `first_sorted_shell` | 形参 |  |
| `jbb` | `middle_sorted_shell` | 形参 |  |
| `jcc` | `last_sorted_shell` | 形参 |  |

### [el1](../../src/lib/spin_angular/two_electron_unchanged_occupations.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `iire` | `check_recoupling` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `ip1` | `first_direct_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip3` | `first_saved_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_direct_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig3` | `last_saved_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `i1` | `two_direct_rank_plus_one` | 局部量 |  |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `i4` | `two_precomputed_rank_plus_one` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `kra1` | `rank_table_position` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `rag` | `subshell_matrix_element` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `pmgg` | `recoupling_by_rank` | 局部量 |  |
| `ragg` | `subshell_element_by_rank` | 局部量 |  |
| `cone` | `exchange_breit_coefficients` | 局部量 |  |

### [el2](../../src/lib/spin_angular/two_electron_pair_transfer.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `iat` | `coupling_allowed` | 局部量 |  |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifaz1` | `orbital_phase_sum` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `jaa` | `first_ordered_shell` | 局部量 |  |
| `jbb` | `last_ordered_shell` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `cond` | `direct_breit_coefficients` | 局部量 |  |

### [el3](../../src/lib/spin_angular/two_electron_single_transfer.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `bra_shell_1` | 形参 |  |
| `jb` | `bra_shell_2` | 形参 |  |
| `jc` | `ket_shell_1` | 形参 |  |
| `jd` | `ket_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |

### [el31](../../src/lib/spin_angular/single_transfer_triple_on_first_shell.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jjja` | `bra_csf_index` | 形参 |  |
| `jjjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jja` | `bra_integral_shell_1` | 形参 |  |
| `jjb` | `bra_integral_shell_2` | 形参 |  |
| `jjc` | `ket_integral_shell_1` | 形参 |  |
| `jjd` | `ket_integral_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `iat` | `coupling_allowed` | 局部量 |  |
| `iia` | `bra_integral_orbital_1` | 局部量 |  |
| `iib` | `bra_integral_orbital_2` | 局部量 |  |
| `iic` | `ket_integral_orbital_1` | 局部量 |  |
| `iid` | `ket_integral_orbital_2` | 局部量 |  |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `jaa` | `first_ordered_shell` | 局部量 |  |
| `jbb` | `last_ordered_shell` | 局部量 |  |
| `jb1` | `phase_end_position` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |

### [el32](../../src/lib/spin_angular/single_transfer_triple_on_second_shell.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jjja` | `bra_csf_index` | 形参 |  |
| `jjjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jja` | `bra_integral_shell_1` | 形参 |  |
| `jjb` | `bra_integral_shell_2` | 形参 |  |
| `jjc` | `ket_integral_shell_1` | 形参 |  |
| `jjd` | `ket_integral_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `iia` | `bra_integral_orbital_1` | 局部量 |  |
| `iib` | `bra_integral_orbital_2` | 局部量 |  |
| `iic` | `ket_integral_orbital_1` | 局部量 |  |
| `iid` | `ket_integral_orbital_2` | 局部量 |  |
| `iatikk` | `unused_rank_status` | 局部量 | 仅声明。 |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `jb1` | `phase_end_position` | 局部量 |  |
| `jaa` | `first_ordered_shell` | 局部量 |  |
| `jbb` | `last_ordered_shell` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |

### [el33](../../src/lib/spin_angular/single_transfer_with_spectator_shell.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jjja` | `bra_csf_index` | 形参 |  |
| `jjjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jc` | `shell_position_c` | 形参 |  |
| `irez` | `orbital_order_case` | 形参 |  |
| `jja` | `bra_integral_shell_1` | 形参 |  |
| `jjb` | `bra_integral_shell_2` | 形参 |  |
| `jjc` | `ket_integral_shell_1` | 形参 |  |
| `jjd` | `ket_integral_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ic` | `orbital_c` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `iia` | `bra_integral_orbital_1` | 局部量 |  |
| `iib` | `bra_integral_orbital_2` | 局部量 |  |
| `iic` | `ket_integral_orbital_1` | 局部量 |  |
| `iid` | `ket_integral_orbital_2` | 局部量 |  |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `i1` | `two_direct_rank_plus_one` | 局部量 |  |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `i4` | `two_precomputed_rank_plus_one` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazp` | `permutation_sign` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `jaa` | `first_ordered_shell` | 局部量 |  |
| `jbb` | `second_ordered_shell` | 局部量 | 先存三壳排序的中间位置，随后用于相位求和时改存相关两个壳层中较后的位置。 |
| `jcc` | `third_ordered_shell` | 局部量 |  |
| `jb1` | `phase_end_before_b` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `kra1` | `rank_table_position` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `l3` | `abs_kappa_3` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `rag` | `subshell_matrix_element` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `pmgg` | `recoupling_by_rank` | 局部量 |  |
| `ragg` | `subshell_element_by_rank` | 局部量 |  |
| `cone` | `exchange_breit_coefficients` | 局部量 |  |

### [el4](../../src/lib/spin_angular/two_electron_three_shell_transfer.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `bra_shell_1` | 形参 |  |
| `jb` | `bra_shell_2` | 形参 |  |
| `jc` | `ket_shell_1` | 形参 |  |
| `jd` | `ket_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |

### [el41](../../src/lib/spin_angular/three_shell_pair_split_or_join.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jjja` | `bra_csf_index` | 形参 |  |
| `jjjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jc` | `shell_position_c` | 形参 |  |
| `irez` | `orbital_order_case` | 形参 |  |
| `jja` | `bra_integral_shell_1` | 形参 |  |
| `jjb` | `bra_integral_shell_2` | 形参 |  |
| `jjc` | `ket_integral_shell_1` | 形参 |  |
| `jjd` | `ket_integral_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ic` | `orbital_c` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `iia` | `bra_integral_orbital_1` | 局部量 |  |
| `iib` | `bra_integral_orbital_2` | 局部量 |  |
| `iic` | `ket_integral_orbital_1` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `iat` | `coupling_allowed` | 局部量 |  |
| `iid` | `ket_integral_orbital_2` | 局部量 |  |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazp` | `permutation_sign` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `jb1` | `phase_end_before_b` | 局部量 |  |
| `jaa` | `first_ordered_shell` | 局部量 |  |
| `jbb` | `second_ordered_shell` | 局部量 | 先存三壳排序的中间位置，随后用于相位求和时改存相关两个壳层中较后的位置。 |
| `jcc` | `third_ordered_shell` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `l3` | `abs_kappa_3` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `cond` | `direct_breit_coefficients` | 局部量 |  |

### [el5](../../src/lib/spin_angular/two_electron_four_shell_transfer.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `bra_shell_1` | 形参 |  |
| `jb` | `bra_shell_2` | 形参 |  |
| `jc` | `ket_shell_1` | 形参 |  |
| `jd` | `ket_shell_2` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |

### [el51](../../src/lib/spin_angular/four_shell_transfer_12_34.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jc` | `shell_position_c` | 形参 |  |
| `jd` | `shell_position_d` | 形参 |  |
| `irez` | `orbital_order_case` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ic` | `orbital_c` | 局部量 |  |
| `id` | `orbital_d` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `i4` | `two_table_rank_plus_one` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazp` | `permutation_sign` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `kra1` | `rank_table_position` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `l3` | `abs_kappa_3` | 局部量 |  |
| `l4` | `abs_kappa_4` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `jb1` | `phase_end_before_b` | 局部量 |  |
| `jd1` | `phase_end_before_d` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `rag` | `subshell_matrix_element` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `pmgg` | `recoupling_by_rank` | 局部量 |  |
| `cond` | `direct_breit_coefficients` | 局部量 |  |
| `cone` | `exchange_breit_coefficients` | 局部量 |  |

### [el52](../../src/lib/spin_angular/four_shell_transfer_13_24.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jc` | `shell_position_c` | 形参 |  |
| `jd` | `shell_position_d` | 形参 |  |
| `irez` | `orbital_order_case` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ic` | `orbital_c` | 局部量 |  |
| `id` | `orbital_d` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `i4` | `two_table_rank_plus_one` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazp` | `permutation_sign` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `kra1` | `rank_table_position` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `l3` | `abs_kappa_3` | 局部量 |  |
| `l4` | `abs_kappa_4` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `jb1` | `phase_end_before_b` | 局部量 |  |
| `jd1` | `phase_end_before_d` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `rag` | `subshell_matrix_element` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `pmgg` | `recoupling_by_rank` | 局部量 |  |
| `cond` | `direct_breit_coefficients` | 局部量 |  |
| `cone` | `exchange_breit_coefficients` | 局部量 |  |

### [el53](../../src/lib/spin_angular/four_shell_transfer_14_23.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `shell_position_a` | 形参 |  |
| `jb` | `shell_position_b` | 形参 |  |
| `jc` | `shell_position_c` | 形参 |  |
| `jd` | `shell_position_d` | 形参 |  |
| `irez` | `orbital_order_case` | 形参 |  |
| `icolbrei` | `interaction_type` | 形参 | 算符选择码。1 为 Coulomb 路径，2 为 Breit 路径。 |
| `ia` | `orbital_a` | 局部量 |  |
| `ib` | `orbital_b` | 局部量 |  |
| `ic` | `orbital_c` | 局部量 |  |
| `id` | `orbital_d` | 局部量 |  |
| `ibrd` | `direct_integral_type` | 局部量 |  |
| `ibre` | `exchange_integral_type` | 局部量 |  |
| `ii` | `work_index` | 局部量 | 根据分支遍历 Breit 工作表或跨壳层费米相位求和区间。 |
| `ip1` | `first_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ip2` | `first_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig1` | `last_coupling_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ig2` | `last_interaction_rank_weight` | 局部量 | 阶数区间的端点，编码为 2k+1。依据参与轨道和耦合次序分组，循环步长为 2。 |
| `ikk` | `two_rank_range_size` | 局部量 |  |
| `i2` | `two_multipole_rank_plus_one` | 局部量 |  |
| `i3` | `two_coupled_rank_plus_one` | 局部量 | 以 2k+1 编码的循环阶数，用于耦合求和。el52 和 el53 也用同一变量遍历直接项。 |
| `i4` | `two_table_rank_plus_one` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `ifazp` | `permutation_sign` | 局部量 |  |
| `ifazfrcs` | `subshell_order_sign` | 局部量 |  |
| `inn` | `phase_orbital` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `kra` | `multipole_rank` | 局部量 |  |
| `kra1` | `rank_table_position` | 局部量 |  |
| `l1` | `abs_kappa_1` | 局部量 |  |
| `l2` | `abs_kappa_2` | 局部量 |  |
| `l3` | `abs_kappa_3` | 局部量 |  |
| `l4` | `abs_kappa_4` | 局部量 |  |
| `j12` | `coupled_rank` | 局部量 |  |
| `jb1` | `phase_end_before_b` | 局部量 |  |
| `jd1` | `phase_end_before_d` | 局部量 |  |
| `nd1` | `first_direct_nu` | 局部量 |  |
| `nd2` | `direct_nu_count` | 局部量 |  |
| `ne1` | `first_exchange_nu` | 局部量 |  |
| `ne2` | `exchange_nu_count` | 局部量 |  |
| `n` | `breit_rank_index` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `nu` | `breit_radial_rank` | 局部量 |  |
| `nup1` | `next_breit_radial_rank` | 局部量 |  |
| `mu` | `breit_component` | 局部量 |  |
| `j` | `two_orbital_j` | 局部量 |  |
| `is` | `integral_orbitals` | 局部量 |  |
| `kaps` | `two_kappas` | 局部量 |  |
| `ks` | `two_j_plus_one` | 局部量 |  |
| `aa` | `coupled_term` | 局部量 |  |
| `ab` | `angular_sum` | 局部量 |  |
| `a1` | `coulomb_angular_factor` | 局部量 |  |
| `bb` | `angular_coefficient` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `qm3` | `q_projection_3` | 局部量 |  |
| `qm4` | `q_projection_4` | 局部量 |  |
| `rag` | `subshell_matrix_element` | 局部量 |  |
| `recc` | `recoupling_factor` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `s` | `breit_factors` | 局部量 |  |
| `pmgg` | `recoupling_by_rank` | 局部量 |  |
| `cond` | `direct_breit_coefficients` | 局部量 |  |
| `cone` | `exchange_breit_coefficients` | 局部量 |  |

### [gg1112](../../src/lib/spin_angular/subshell_product_wa_a.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik1` | `bra_shell_codes_1` | 形参 |  |
| `ik2` | `bra_shell_codes_2` | 形参 |  |
| `id1` | `ket_shell_codes_1` | 形参 |  |
| `id2` | `ket_shell_codes_2` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `qm4` | `q_projection_4` | 形参 |  |
| `bk1` | `bra_shell_values_1` | 形参 |  |
| `bk2` | `bra_shell_values_2` | 形参 |  |
| `bd1` | `ket_shell_values_1` | 形参 |  |
| `bd2` | `ket_shell_values_2` | 形参 |  |
| `ww` | `subshell_matrix_element` | 形参 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `iqmm3` | `occupation_change_3` | 局部量 |  |
| `iqmm4` | `occupation_change_4` | 局部量 |  |
| `iqmm23` | `first_shell_occupation_change` | 局部量 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `a1` | `cg_factor` | 局部量 |  |
| `s` | `reduced_cfp` | 局部量 |  |
| `bk` | `three_operator_rank` | 局部量 |  |
| `wa` | `wa_matrix_element` | 局部量 |  |

### [gg1122](../../src/lib/spin_angular/subshell_product_w_w.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `first_w_rank` | 形参 |  |
| `k2` | `second_w_rank` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `qm4` | `q_projection_4` | 形参 |  |
| `aa` | `subshell_matrix_element` | 形参 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `iqmm3` | `occupation_change_3` | 局部量 |  |
| `iqmm4` | `occupation_change_4` | 局部量 |  |
| `iqmm12` | `first_shell_occupation_change` | 局部量 |  |
| `iqmm34` | `last_shell_occupation_change` | 局部量 |  |
| `a1` | `first_w_matrix_element` | 局部量 |  |
| `w` | `second_w_matrix_element` | 局部量 |  |

### [gg12](../../src/lib/spin_angular/subshell_product_a_a.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ik1` | `bra_shell_codes_1` | 形参 |  |
| `ik2` | `bra_shell_codes_2` | 形参 |  |
| `id1` | `ket_shell_codes_1` | 形参 |  |
| `id2` | `ket_shell_codes_2` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `bk1` | `bra_shell_values_1` | 形参 |  |
| `bk2` | `bra_shell_values_2` | 形参 |  |
| `bd1` | `ket_shell_values_1` | 形参 |  |
| `bd2` | `ket_shell_values_2` | 形参 |  |
| `ww` | `subshell_matrix_element` | 形参 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `a1` | `cg_product` | 局部量 |  |
| `c` | `second_reduced_cfp` | 局部量 |  |
| `c1` | `second_cg_factor` | 局部量 |  |
| `s` | `first_reduced_cfp` | 局部量 |  |

### [gg1222](../../src/lib/spin_angular/subshell_product_a_aw.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik1` | `bra_shell_codes_1` | 形参 |  |
| `ik2` | `bra_shell_codes_2` | 形参 |  |
| `id1` | `ket_shell_codes_1` | 形参 |  |
| `id2` | `ket_shell_codes_2` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `qm4` | `q_projection_4` | 形参 |  |
| `bk1` | `bra_shell_values_1` | 形参 |  |
| `bk2` | `bra_shell_values_2` | 形参 |  |
| `bd1` | `ket_shell_values_1` | 形参 |  |
| `bd2` | `ket_shell_values_2` | 形参 |  |
| `ww` | `subshell_matrix_element` | 形参 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `iqmm3` | `occupation_change_3` | 局部量 |  |
| `iqmm4` | `occupation_change_4` | 局部量 |  |
| `iqmm34` | `last_shell_occupation_change` | 局部量 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `a1` | `cg_factor` | 局部量 |  |
| `s` | `reduced_cfp` | 局部量 |  |
| `bk` | `three_operator_rank` | 局部量 |  |
| `aw` | `aw_matrix_element` | 局部量 |  |

### [gg1233](../../src/lib/spin_angular/subshell_product_a_a_w.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik1` | `bra_shell_codes_1` | 形参 |  |
| `ik2` | `bra_shell_codes_2` | 形参 |  |
| `ik3` | `bra_shell_codes_3` | 形参 |  |
| `id1` | `ket_shell_codes_1` | 形参 |  |
| `id2` | `ket_shell_codes_2` | 形参 |  |
| `id3` | `ket_shell_codes_3` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `qm4` | `q_projection_4` | 形参 |  |
| `bk1` | `bra_shell_values_1` | 形参 |  |
| `bk2` | `bra_shell_values_2` | 形参 |  |
| `bk3` | `bra_shell_values_3` | 形参 |  |
| `bd1` | `ket_shell_values_1` | 形参 |  |
| `bd2` | `ket_shell_values_2` | 形参 |  |
| `bd3` | `ket_shell_values_3` | 形参 |  |
| `ww` | `subshell_matrix_element` | 形参 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `iqmm3` | `occupation_change_3` | 局部量 |  |
| `iqmm4` | `occupation_change_4` | 局部量 |  |
| `iqmm34` | `last_shell_occupation_change` | 局部量 |  |
| `a1` | `cg_product` | 局部量 |  |
| `c` | `second_reduced_cfp` | 局部量 |  |
| `c1` | `next_cg_factor` | 局部量 |  |
| `s` | `first_reduced_cfp` | 局部量 |  |
| `w` | `last_w_matrix_element` | 局部量 |  |

### [gg1234](../../src/lib/spin_angular/subshell_product_a_a_a_a.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ik1` | `bra_shell_codes_1` | 形参 |  |
| `ik2` | `bra_shell_codes_2` | 形参 |  |
| `ik3` | `bra_shell_codes_3` | 形参 |  |
| `ik4` | `bra_shell_codes_4` | 形参 |  |
| `id1` | `ket_shell_codes_1` | 形参 |  |
| `id2` | `ket_shell_codes_2` | 形参 |  |
| `id3` | `ket_shell_codes_3` | 形参 |  |
| `id4` | `ket_shell_codes_4` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `qm4` | `q_projection_4` | 形参 |  |
| `bk1` | `bra_shell_values_1` | 形参 |  |
| `bk2` | `bra_shell_values_2` | 形参 |  |
| `bk3` | `bra_shell_values_3` | 形参 |  |
| `bk4` | `bra_shell_values_4` | 形参 |  |
| `bd1` | `ket_shell_values_1` | 形参 |  |
| `bd2` | `ket_shell_values_2` | 形参 |  |
| `bd3` | `ket_shell_values_3` | 形参 |  |
| `bd4` | `ket_shell_values_4` | 形参 |  |
| `ww` | `subshell_matrix_element` | 形参 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `iqmm3` | `occupation_change_3` | 局部量 |  |
| `iqmm4` | `occupation_change_4` | 局部量 |  |
| `a1` | `cg_product` | 局部量 |  |
| `c` | `second_reduced_cfp` | 局部量 |  |
| `c1` | `next_cg_factor` | 局部量 |  |
| `s` | `first_reduced_cfp` | 局部量 |  |
| `v` | `third_reduced_cfp` | 局部量 |  |
| `z` | `fourth_reduced_cfp` | 局部量 |  |

### [gracah1](../../src/lib/spin_angular/racah_coefficients_scott_fischer.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `two_j_1` | 形参 |  |
| `j` | `two_j_2` | 形参 |  |
| `k` | `two_j_3` | 形参 |  |
| `l` | `two_j_4` | 形参 |  |
| `m` | `two_j_5` | 形参 |  |
| `n` | `two_j_6` | 形参 |  |
| `rac` | `racah_coefficient` | 形参 |  |
| `j1` | `triangle_sum_1` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j2` | `triangle_sum_2` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j3` | `triangle_sum_3` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j4` | `triangle_sum_4` | 局部量 | 三角组的角动量和。先保存二倍值供选择规则检查，再除以 2。 |
| `j5` | `quadrangle_sum_1` | 局部量 |  |
| `j6` | `quadrangle_sum_2` | 局部量 |  |
| `j7` | `quadrangle_sum_3` | 局部量 |  |
| `numin` | `sum_start` | 局部量 |  |
| `numax` | `sum_end` | 局部量 |  |
| `icount` | `terms_evaluated` | 局部量 |  |
| `kk` | `term_iteration` | 局部量 |  |
| `ki` | `descending_term_index` | 局部量 |  |

### [itjj](../../src/lib/spin_angular/jj_intermediate_terms_same_parity.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `kg` | `tensor_rank` | 形参 |  |
| `iq` | `occupation_change` | 形参 |  |
| `kg1` | `two_tensor_rank` | 形参 |  |
| `itp` | `first_intermediate_term` | 形参 |  |
| `itg` | `last_intermediate_term` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `ibt` | `intermediate_shell_codes` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `bt` | `intermediate_shell_values` | 形参 |  |
| `itk` | `bra_term_lookup_index` | 局部量 |  |
| `itd` | `ket_term_lookup_index` | 局部量 |  |
| `itp1` | `first_bra_range_term` | 局部量 |  |
| `itg1` | `last_bra_range_term` | 局部量 |  |
| `itjj` | `coupling_allowed` | 函数结果 |  |

### [itjj2](../../src/lib/spin_angular/jj_intermediate_terms_ket_parity.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `kg1` | `two_tensor_rank` | 形参 |  |
| `iq` | `occupation_change` | 形参 |  |
| `itp` | `first_intermediate_term` | 形参 |  |
| `itg` | `last_intermediate_term` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `ibt` | `intermediate_shell_codes` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `bt` | `intermediate_shell_values` | 形参 |  |
| `itk` | `bra_term_lookup_index` | 局部量 |  |
| `itd` | `ket_term_lookup_index` | 局部量 |  |
| `itp1` | `first_bra_range_term` | 局部量 |  |
| `itg1` | `last_bra_range_term` | 局部量 |  |
| `itjj2` | `coupling_allowed` | 函数结果 |  |

### [itjj3](../../src/lib/spin_angular/jj_intermediate_terms_bra_parity.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `kg1` | `two_tensor_rank` | 形参 |  |
| `iq` | `occupation_change` | 形参 |  |
| `itp` | `first_intermediate_term` | 形参 |  |
| `itg` | `last_intermediate_term` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `ibt` | `intermediate_shell_codes` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `bt` | `intermediate_shell_values` | 形参 |  |
| `itk` | `bra_term_lookup_index` | 局部量 |  |
| `itd` | `ket_term_lookup_index` | 局部量 |  |
| `itp1` | `first_bra_range_term` | 局部量 |  |
| `itg1` | `last_bra_range_term` | 局部量 |  |
| `itjj3` | `coupling_allowed` | 函数结果 |  |

### [itrexg](../../src/lib/spin_angular/shared_tensor_rank_range.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i1` | `two_j_1` | 形参 |  |
| `i2` | `two_j_2` | 形参 |  |
| `i3` | `two_j_3` | 形参 |  |
| `i4` | `two_j_4` | 形参 |  |
| `k` | `two_rank_range_size` | 形参 | 闭区间 max_two_rank - min_two_rank + 1，调用方使用步长 2 选取实际允许的阶数。 |
| `j` | `lower_two_rank` | 局部量 |  |
| `itrexg` | `min_two_rank` | 函数结果 |  |

### [ittk](../../src/lib/spin_angular/angular_triangle_and_parity.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `first_angular_label` | 形参 | 调用方传入统一编码。Wigner 路径传 2j，coulom 路径传整数 l 和 k 以同时检查奇偶。 |
| `j` | `second_angular_label` | 形参 |  |
| `k` | `third_angular_label` | 形参 |  |
| `ittk` | `triangle_allowed` | 函数结果 |  |

### [ixjtik](../../src/lib/spin_angular/wigner_6j_triangle_selection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `two_j_1` | 形参 |  |
| `j` | `two_j_2` | 形参 |  |
| `k` | `two_j_3` | 形参 |  |
| `l` | `two_j_4` | 形参 |  |
| `m` | `two_j_5` | 形参 |  |
| `n` | `two_j_6` | 形参 |  |
| `ixjtik` | `triangles_allowed` | 函数结果 |  |

### [izas1](../../src/lib/spin_angular/quasispin_projection_selection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ib` | `two_ket_q` | 形参 |  |
| `ik` | `two_bra_q` | 形参 |  |
| `qb` | `ket_q_projection` | 形参 |  |
| `qk` | `bra_q_projection` | 形参 |  |
| `iqb` | `two_abs_ket_projection` | 局部量 |  |
| `iqk` | `two_abs_bra_projection` | 局部量 |  |
| `izas1` | `projections_allowed` | 函数结果 |  |

### [jfaze](../../src/lib/spin_angular/subshell_permutation_sign.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i1` | `first_shell_position` | 形参 |  |
| `i2` | `second_shell_position` | 形参 |  |
| `i3` | `third_shell_position` | 形参 |  |
| `i4` | `fourth_shell_position` | 形参 |  |
| `jfaze` | `permutation_sign` | 函数结果 |  |

### [jthn](../../src/lib/spin_angular/packed_integer_digit.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `radix` | 形参 |  |
| `n` | `digit_position` | 形参 |  |
| `k` | `packed_value` | 形参 |  |
| `jthn` | `digit_value` | 函数结果 |  |

### [mes](../../src/lib/spin_angular/subshell_tensor_diagnostics.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `diagnostic_code` | 形参 |  |
| `j` | `message_index` | 局部量 |  |
| `string5` | `routine_names` | DATA 初始化 |  |

### [nine](../../src/lib/spin_angular/wigner_9j.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `two_j_11` | 形参 |  |
| `j2` | `two_j_12` | 形参 |  |
| `j3` | `two_j_13` | 形参 |  |
| `l1` | `two_j_21` | 形参 |  |
| `l2` | `two_j_22` | 形参 |  |
| `l3` | `two_j_23` | 形参 |  |
| `k1` | `two_j_31` | 形参 |  |
| `k2` | `two_j_32` | 形参 |  |
| `k3` | `two_j_33` | 形参 |  |
| `i` | `check_triangles_only` | 形参 |  |
| `in` | `evaluation_status` | 形参 | 只做三角条件检查时表示通过与否。求值路径置 1，数值由 nine_j_value 返回。 |
| `aa` | `nine_j_value` | 形参 |  |
| `n1` | `triangle_bound_1` | 局部量 |  |
| `n2` | `triangle_bound_2` | 局部量 |  |
| `n3` | `triangle_bound_3` | 局部量 |  |
| `n4` | `triangle_bound_4` | 局部量 |  |
| `n5` | `triangle_bound_5` | 局部量 |  |
| `n6` | `triangle_bound_6` | 局部量 |  |
| `max_` | `min_two_intermediate_j` | 局部量 |  |
| `min_` | `max_two_intermediate_j` | 局部量 |  |
| `ix` | `two_intermediate_j` | 局部量 |  |
| `s1` | `first_six_j` | 局部量 |  |
| `s2` | `second_six_j` | 局部量 |  |
| `s3` | `third_six_j` | 局部量 |  |
| `x` | `weighted_six_j_product` | 局部量 |  |

### [nine0](../../src/lib/spin_angular/wigner_9j_with_zero.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `two_j_11` | 形参 |  |
| `j2` | `two_j_12` | 形参 |  |
| `j3` | `two_j_13` | 形参 |  |
| `l1` | `two_j_21` | 形参 |  |
| `l2` | `two_j_22` | 形参 |  |
| `l3` | `two_j_23` | 形参 |  |
| `k1` | `two_j_31` | 形参 |  |
| `k2` | `two_j_32` | 形参 |  |
| `k3` | `two_j_33` | 形参 |  |
| `aa` | `nine_j_value` | 形参 |  |
| `ifa` | `twice_phase_exponent` | 局部量 |  |
| `a` | `six_j_value` | 局部量 |  |
| `b` | `normalization_squared` | 局部量 |  |

### [nmtejj](../../src/lib/spin_angular/jj_subshell_term_lookup.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i2q` | `two_q` | 形参 |  |
| `i2j` | `two_term_j` | 形参 |  |
| `j` | `two_orbital_j` | 形参 |  |
| `nk` | `ket_occupation` | 形参 |  |
| `nd` | `bra_occupation` | 形参 |  |
| `i2v` | `two_seniority` | 局部量 |  |
| `inn` | `packed_term_code` | 局部量 |  |
| `il` | `high_j_range_index` | 局部量 |  |
| `jp` | `first_term_index` | 局部量 |  |
| `jg` | `last_term_index` | 局部量 |  |
| `jj` | `table_index` | 局部量 |  |
| `lp` | `low_j_first_terms` | DATA 初始化 |  |
| `lg` | `low_j_last_terms` | DATA 初始化 |  |
| `lp3` | `high_j_first_terms` | DATA 初始化 |  |
| `lg3` | `high_j_last_terms` | DATA 初始化 |  |
| `nmtejj` | `term_index` | 函数结果 |  |

### [oneparticlejj](../../src/lib/spin_angular/one_electron_tensor_coefficients.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ka` | `tensor_rank` | 形参 |  |
| `iopar` | `operator_parity` | 形参 |  |
| `ja` | `bra_csf_index` | 形参 |  |
| `jb` | `ket_csf_index` | 形参 |  |
| `ia1` | `bra_orbital_index` | 形参 |  |
| `ia2` | `ket_orbital_index` | 形参 |  |
| `vshell` | `shell_coefficients` | 形参 | 占据相同时用完整轨道号索引。单电子转移时仅第一项返回系数，两个轨道号单独输出。 |
| `kk` | `two_rank_plus_one` | 局部量 |  |
| `ij` | `orbital_index` | 局部量 |  |
| `idq` | `total_occupation_difference` | 局部量 | 两个 CSF 子壳层占据差的绝对值之和。单电子转移时为 2。 |
| `ja1` | `bra_shell_position` | 局部量 |  |
| `ja2` | `ket_shell_position` | 局部量 |  |
| `jw` | `peel_position` | 局部量 |  |
| `ndq` | `occupation_difference` | 局部量 |  |
| `ns` | `active_shell_count` | 局部量 |  |
| `ish` | `diagonal_orbital` | 局部量 |  |
| `i` | `shell_position` | 局部量 |  |
| `ii` | `shift_step` | 局部量 |  |
| `i1` | `first_orbital` | 局部量 |  |
| `im` | `shift_count` | 局部量 |  |
| `ks1` | `two_bra_j_plus_one` | 局部量 |  |
| `ks2` | `two_ket_j_plus_one` | 局部量 |  |
| `nx` | `orbital_difference` | 局部量 |  |
| `npeelm` | `peel_coupling_count` | 局部量 |  |
| `is` | `unused_orbital_labels` | 局部量 | 仅声明。 |
| `ks` | `unused_orbital_j_weights` | 局部量 | 仅声明。 |
| `tcoeff` | `angular_coefficient` | 局部量 |  |

### [oneparticlejj1](../../src/lib/spin_angular/one_electron_tensor_same_subshell.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `ka` | `tensor_rank` | 形参 |  |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `bra_shell_position` | 形参 |  |
| `jb` | `ket_shell_position` | 形参 |  |
| `coeff` | `angular_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `rec` | `recoupling_factor` | 局部量 |  |
| `wj` | `w_matrix_element` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |

### [oneparticlejj2](../../src/lib/spin_angular/one_electron_tensor_between_subshells.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `ka` | `tensor_rank` | 形参 |  |
| `ja` | `bra_shell_position` | 形参 |  |
| `jb` | `ket_shell_position` | 形参 |  |
| `coeff` | `angular_coefficient` | 形参 |  |
| `i` | `shell_position` | 局部量 |  |
| `ia` | `first_shell_position` | 局部量 |  |
| `ib` | `last_shell_position` | 局部量 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `ij` | `orbital_index` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `iqmm1` | `occupation_change_1` | 局部量 |  |
| `iqmm2` | `occupation_change_2` | 局部量 |  |
| `jibks1` | `unused_shell_workspace` | 局部量 | 仅声明。 |
| `jib` | `last_phase_shell` | 局部量 |  |
| `ks1` | `two_first_orbital_j` | 局部量 |  |
| `ks2` | `two_last_orbital_j` | 局部量 |  |
| `rec` | `recoupling_factor` | 局部量 |  |
| `a1` | `subshell_matrix_element` | 局部量 |  |
| `a2` | `first_cg_factor` | 局部量 |  |
| `a3` | `second_cg_factor` | 局部量 |  |
| `s1` | `first_reduced_cfp` | 局部量 |  |
| `s2` | `second_reduced_cfp` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |

### [onescalar](../../src/lib/spin_angular/one_electron_scalar_coefficients.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja` | `bra_csf_index` | 形参 |  |
| `jb` | `ket_csf_index` | 形参 |  |
| `ia1` | `bra_orbital_index` | 形参 |  |
| `ia2` | `ket_orbital_index` | 形参 |  |
| `vshell` | `shell_coefficients` | 形参 | 占据相同时用完整轨道号索引。单电子转移时仅第一项返回系数，两个轨道号单独输出。 |
| `kk` | `two_rank_plus_one` | 局部量 |  |
| `iopar` | `operator_parity` | 局部量 |  |
| `ij` | `orbital_index` | 局部量 |  |
| `idq` | `total_occupation_difference` | 局部量 | 两个 CSF 子壳层占据差的绝对值之和。单电子转移时为 2。 |
| `ja1` | `bra_shell_position` | 局部量 |  |
| `ja2` | `ket_shell_position` | 局部量 |  |
| `ndq` | `occupation_difference` | 局部量 |  |
| `ns` | `active_shell_count` | 局部量 |  |
| `ish` | `diagonal_orbital` | 局部量 |  |
| `i` | `shell_position` | 局部量 |  |
| `ii` | `shift_step` | 局部量 |  |
| `i1` | `first_orbital` | 局部量 |  |
| `im` | `shift_count` | 局部量 |  |
| `jw` | `peel_position` | 局部量 |  |
| `ks1` | `two_bra_j_plus_one` | 局部量 |  |
| `ks2` | `two_ket_j_plus_one` | 局部量 |  |
| `nx` | `orbital_difference` | 局部量 |  |
| `npeelm` | `peel_coupling_count` | 局部量 |  |
| `is` | `unused_orbital_labels` | 局部量 | 仅声明。 |
| `ks` | `unused_orbital_j_weights` | 局部量 | 仅声明。 |
| `tcoeff` | `angular_coefficient` | 局部量 |  |

### [onescalar1](../../src/lib/spin_angular/one_electron_scalar_same_subshell.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `bra_shell_position` | 形参 |  |
| `jb` | `ket_shell_position` | 形参 |  |
| `coeff` | `angular_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `wj` | `w_matrix_element` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `recoupl` | `recoupling_factor` | 局部量 |  |

### [onescalar2](../../src/lib/spin_angular/one_electron_scalar_between_subshells.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jja` | `bra_csf_index` | 形参 |  |
| `jjb` | `ket_csf_index` | 形参 |  |
| `ja` | `bra_shell_position` | 形参 |  |
| `jb` | `ket_shell_position` | 形参 |  |
| `coeff` | `angular_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 局部量 |  |
| `jaa` | `first_shell_position` | 局部量 |  |
| `jbb` | `last_shell_position` | 局部量 |  |
| `nn` | `intervening_electron_count` | 局部量 |  |
| `ib1` | `last_phase_shell` | 局部量 |  |
| `ii` | `shell_position` | 局部量 |  |
| `inn` | `orbital_index` | 局部量 |  |
| `ifaz` | `phase_sum` | 局部量 |  |
| `qm1` | `q_projection_1` | 局部量 |  |
| `qm2` | `q_projection_2` | 局部量 |  |
| `rec` | `recoupling_factor` | 局部量 |  |
| `ww` | `subshell_matrix_element` | 局部量 |  |

### [perko1](../../src/lib/spin_angular/jj_subshell_state_packing.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja` | `shell_position` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `ij` | `orbital_index` | 局部量 |  |

### [perko2](../../src/lib/spin_angular/jj_active_subshell_states.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ja3` | `third_shell_position` | 形参 |  |
| `ja4` | `fourth_shell_position` | 形参 |  |
| `i` | `active_shell_count` | 形参 |  |

### [rec3](../../src/lib/spin_angular/permuted_three_subshell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ja3` | `third_shell_position` | 形参 |  |
| `k1` | `two_rank_1` | 形参 |  |
| `k2` | `two_rank_2` | 形参 |  |
| `k3` | `two_rank_3` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ifaz` | `phase_sum` | 局部量 |  |

### [reco](../../src/lib/spin_angular/scalar_interaction_recoupling_selection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ja3` | `third_shell_position` | 形参 |  |
| `ja4` | `fourth_shell_position` | 形参 |  |
| `ka` | `interaction_case` | 形参 | 决定检查时跳过哪些活跃子壳层及中间耦合位置的分支编号。 |
| `iat` | `coupling_allowed` | 形参 |  |
| `i` | `shell_position` | 局部量 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `ij` | `orbital_index` | 局部量 |  |
| `j` | `shell_code_slot` | 局部量 |  |
| `ia1` | `first_skipped_coupling` | 局部量 |  |
| `ia2` | `coupling_skip_end` | 局部量 | 跳过的中间耦合区间上界，采用开区间端点。 |

### [reco2](../../src/lib/spin_angular/two_subshell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ia1` | `two_bra_shell_j_1` | 局部量 |  |
| `ia2` | `two_bra_shell_j_2` | 局部量 |  |
| `ib1` | `two_ket_shell_j_1` | 局部量 |  |
| `ib2` | `two_ket_shell_j_2` | 局部量 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `iskr` | `coupling_step_count` | 局部量 |  |
| `s` | `first_bra_j_weight` | 局部量 |  |
| `ss` | `bra_j_weight_product` | 局部量 |  |
| `re` | `diagram_factor` | 局部量 |  |

### [reco3](../../src/lib/spin_angular/three_subshell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ja3` | `third_shell_position` | 形参 |  |
| `k1` | `two_rank_1` | 形参 |  |
| `k2` | `two_rank_2` | 形参 |  |
| `ka` | `two_tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ia3` | `two_bra_shell_j_3` | 局部量 |  |
| `ib3` | `two_ket_shell_j_3` | 局部量 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `ij3` | `orbital_index_3` | 局部量 |  |
| `iskr` | `coupling_step_count` | 局部量 |  |
| `s` | `bra_j_weight_product` | 局部量 |  |
| `s1` | `bra_j_weight_1` | 局部量 |  |
| `s2` | `bra_j_weight_2` | 局部量 |  |
| `s3` | `bra_j_weight_3` | 局部量 |  |
| `re` | `diagram_factor` | 局部量 |  |

### [reco4](../../src/lib/spin_angular/four_subshell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ja3` | `third_shell_position` | 形参 |  |
| `ja4` | `fourth_shell_position` | 形参 |  |
| `k1` | `two_rank_1` | 形参 |  |
| `k2` | `two_rank_2` | 形参 |  |
| `k3` | `two_rank_3` | 形参 |  |
| `k4` | `two_second_pair_rank` | 形参 |  |
| `ka` | `two_first_pair_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ia4` | `two_bra_shell_j_4` | 局部量 |  |
| `ib4` | `two_ket_shell_j_4` | 局部量 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `ij3` | `orbital_index_3` | 局部量 |  |
| `ij4` | `orbital_index_4` | 局部量 |  |
| `iskr` | `coupling_step_count` | 局部量 |  |
| `s` | `bra_j_weight_product` | 局部量 |  |
| `s1` | `bra_j_weight_1` | 局部量 |  |
| `s2` | `bra_j_weight_2` | 局部量 |  |
| `s3` | `bra_j_weight_3` | 局部量 |  |
| `s4` | `bra_j_weight_4` | 局部量 |  |
| `re` | `diagram_factor` | 局部量 |  |

### [recoonescalar](../../src/lib/spin_angular/one_electron_scalar_recoupling_selection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ja3` | `third_shell_position` | 形参 |  |
| `ja4` | `fourth_shell_position` | 形参 |  |
| `ka` | `interaction_case` | 形参 | 决定检查时跳过哪些活跃子壳层及中间耦合位置的分支编号。 |
| `iat` | `coupling_allowed` | 形参 |  |
| `i` | `shell_position` | 局部量 |  |
| `ia1` | `first_skipped_coupling` | 局部量 |  |
| `ia2` | `coupling_skip_end` | 局部量 | 跳过的中间耦合区间上界，采用开区间端点。 |
| `ij` | `orbital_index` | 局部量 |  |
| `ij1` | `orbital_index_1` | 局部量 |  |
| `ij2` | `orbital_index_2` | 局部量 |  |
| `j` | `shell_code_slot` | 局部量 |  |
| `npeelgg` | `coupled_shell_count` | 局部量 |  |

### [recop00](../../src/lib/spin_angular/one_electron_tensor_recoupling_selection.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `ka` | `tensor_rank` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `i` | `shell_position` | 局部量 |  |
| `iskr` | `coupling_step_count` | 局部量 |  |
| `ij1` | `orbital_index` | 局部量 |  |
| `j` | `prefix_shell_position` | 局部量 |  |
| `ji` | `suffix_step` | 局部量 |  |
| `jj` | `prefix_coupling_position` | 局部量 |  |
| `npeelgg` | `coupled_shell_count` | 局部量 |  |
| `kk` | `suffix_coupling_position` | 局部量 |  |
| `lk` | `two_bra_chain_j` | 局部量 |  |
| `ld` | `two_ket_chain_j` | 局部量 |  |

### [recop1](../../src/lib/spin_angular/one_electron_same_subshell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `ja1` | `first_shell_position` | 形参 |  |
| `ka` | `tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ij1` | `orbital_index` | 局部量 |  |
| `iskr` | `coupling_step_count` | 局部量 |  |
| `npeelgg` | `coupled_shell_count` | 局部量 |  |
| `s` | `bra_subshell_j_weight` | 局部量 |  |
| `re` | `diagram_factor` | 局部量 |  |

### [recop2](../../src/lib/spin_angular/one_electron_two_subshell_recoupling.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ns` | `active_shell_count` | 形参 |  |
| `ja1` | `first_shell_position` | 形参 |  |
| `ja2` | `second_shell_position` | 形参 |  |
| `k1` | `two_rank_1` | 形参 |  |
| `k2` | `two_rank_2` | 形参 |  |
| `ka` | `tensor_rank` | 形参 |  |
| `ire` | `evaluate_coefficient` | 形参 |  |
| `iat` | `coupling_allowed` | 形参 |  |
| `recc` | `recoupling_coefficient` | 形参 |  |
| `ij` | `orbital_index` | 局部量 |  |
| `iskr` | `coupling_step_count` | 局部量 |  |
| `npeelgg` | `coupled_shell_count` | 局部量 |  |
| `s` | `first_bra_j_weight` | 局部量 |  |
| `ss` | `bra_j_weight_product` | 局部量 |  |
| `re` | `diagram_factor` | 局部量 |  |

### [rkco_gg](../../src/lib/spin_angular/two_electron_angular_coefficients.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `coefficient_consumer` | `coefficient_consumer` | 过程形参 | 处理闭壳层贡献的过程形参。 |
| `ja` | `bra_csf_index` | 形参 |  |
| `jb` | `ket_csf_index` | 形参 |  |
| `incor` | `include_closed_shell_terms` | 形参 | 小于 1 时跳过文件末端闭壳层回调。开壳层与闭壳层交叉项仍由前段路径处理。 |
| `icolbrei` | `interaction_type` | 形参 |  |
| `i` | `shell_position` | 局部量 |  |
| `ib1` | `bra_orbital` | 局部量 |  |
| `ib2` | `ket_orbital` | 局部量 |  |
| `idq` | `total_occupation_difference` | 局部量 |  |
| `idqg` | `occupation_pattern_code` | 局部量 | 每个相差 1 电子的子壳层加 1，每个相差 2 电子的子壳层加 20，用于分派不同占据差模式。 |
| `ii` | `shift_step` | 局部量 |  |
| `i1` | `first_orbital` | 局部量 |  |
| `ij` | `orbital_index` | 局部量 |  |
| `ijw` | `core_orbital` | 局部量 |  |
| `im` | `shift_count` | 局部量 |  |
| `ipca` | `orbital_list_selector` | 局部量 |  |
| `it1` | `first_bra_orbital` | 局部量 |  |
| `it2` | `first_ket_orbital` | 局部量 |  |
| `ja1` | `bra_shell_1` | 局部量 | 双电子贡献的左侧第一个子壳层位置。闭壳层回调按 orbital_list_selector 解释其位置。 |
| `jb1` | `bra_shell_2` | 局部量 |  |
| `ja2` | `ket_shell_1` | 局部量 |  |
| `jb2` | `ket_shell_2` | 局部量 |  |
| `j` | `peel_orbital` | 局部量 |  |
| `jw` | `peel_position` | 局部量 |  |
| `jt1` | `saved_bra_shell_1` | 局部量 |  |
| `jt2` | `saved_ket_shell_1` | 局部量 |  |
| `jt3` | `swap_shell_position` | 局部量 |  |
| `ndq` | `occupation_difference` | 局部量 |  |
| `klast` | `last_outer_position` | 局部量 |  |
| `kw` | `partner_shell_position` | 局部量 |  |
| `kwa` | `outer_peel_position` | 局部量 |  |
| `k1` | `partner_orbital` | 局部量 |  |
| `kw1` | `second_shell_position` | 局部量 |  |
| `kw2` | `first_shell_position` | 局部量 |  |
| `npeelm` | `peel_coupling_count` | 局部量 |  |

### [rmeajj](../../src/lib/spin_angular/jj_reduced_cfp.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ll` | `two_orbital_j` | 形参 |  |
| `it` | `bra_term_index` | 形参 |  |
| `lq` | `two_bra_q` | 形参 |  |
| `j` | `two_bra_term_j` | 形参 |  |
| `its` | `ket_term_index` | 形参 |  |
| `lqs` | `two_ket_q` | 形参 |  |
| `j1s` | `two_ket_term_j` | 形参 |  |
| `coef` | `reduced_cfp` | 形参 |  |
| `i1` | `cfp_row` | 局部量 |  |
| `i2` | `cfp_column` | 局部量 |  |
| `ie1` | `twice_phase_exponent` | 局部量 |  |
| `jt` | `lower_term_index` | 局部量 |  |
| `jts` | `upper_term_index` | 局部量 |  |
| `ids3` | `cfp_numerators_two_j_3` | DATA 初始化 |  |
| `idv3` | `cfp_denominators_two_j_3` | DATA 初始化 |  |
| `ids5` | `cfp_numerators_two_j_5` | DATA 初始化 |  |
| `idv5` | `cfp_denominators_two_j_5` | DATA 初始化 |  |
| `ids7` | `cfp_numerators_two_j_7` | DATA 初始化 |  |
| `idv7` | `cfp_denominators_two_j_7` | DATA 初始化 |  |

### [rmeajj11](../../src/lib/spin_angular/jj_reduced_cfp_two_electron_terms.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `ll` | `two_orbital_j` | 形参 |  |
| `j1` | `two_electron_term_index` | 形参 |  |
| `j2` | `one_electron_term_index` | 形参 |  |
| `s` | `reduced_cfp` | 形参 |  |
| `lq` | `two_bra_q` | 局部量 |  |
| `lv` | `two_bra_seniority` | 局部量 |  |
| `lqs` | `two_ket_q` | 局部量 |  |
| `lvs` | `two_ket_seniority` | 局部量 |  |
| `j` | `two_bra_term_j` | 局部量 |  |
| `j1s` | `two_ket_term_j` | 局部量 |  |
| `n` | `bra_occupation` | 局部量 |  |
| `inn` | `phase_exponent` | 局部量 |  |
| `a1` | `normalization_squared` | 局部量 |  |
| `a4` | `cg_coefficient` | 局部量 |  |
| `q` | `bra_q` | 局部量 |  |
| `qq` | `operator_q_projection` | 局部量 |  |
| `qs` | `ket_q` | 局部量 |  |
| `qm` | `bra_q_projection` | 局部量 |  |
| `qms` | `ket_q_projection` | 局部量 |  |

### [rmeajj9](../../src/lib/spin_angular/jj_reduced_cfp_j_nine_halves.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `it` | `bra_term_index` | 形参 |  |
| `lq` | `two_bra_q` | 形参 | 接口携带的左侧 2Q。当前表按项索引，函数体不读取该参数。 |
| `j` | `two_bra_term_j` | 形参 | 接口携带的左侧子壳层项 2J。当前表按项索引，函数体不读取该参数。 |
| `its` | `ket_term_index` | 形参 |  |
| `lqs` | `two_ket_q` | 形参 | 接口携带的右侧 2Q。当前表按项索引，函数体不读取该参数。 |
| `j1s` | `two_ket_term_j` | 形参 | 接口携带的右侧子壳层项 2J。当前表按项索引，函数体不读取该参数。 |
| `coef` | `reduced_cfp` | 形参 |  |
| `jt` | `cfp_row` | 局部量 |  |
| `jts` | `cfp_column` | 局部量 |  |
| `ids3` | `cfp_numerators_col_3` | DATA 初始化 |  |
| `ids4` | `cfp_numerators_col_4` | DATA 初始化 |  |
| `ids5` | `cfp_numerators_col_5` | DATA 初始化 |  |
| `ids6` | `cfp_numerators_col_6` | DATA 初始化 |  |
| `ids7` | `cfp_numerators_col_7` | DATA 初始化 |  |
| `idv7` | `cfp_denominators_col_7` | DATA 初始化 |  |
| `ids8` | `cfp_numerators_col_8` | DATA 初始化 |  |
| `idv8` | `cfp_denominators_col_8` | DATA 初始化 |  |
| `ids9` | `cfp_numerators_col_9` | DATA 初始化 |  |
| `ids10` | `cfp_numerators_col_10` | DATA 初始化 |  |
| `ids11` | `cfp_numerators_col_11` | DATA 初始化 |  |
| `idv11` | `cfp_denominators_col_11` | DATA 初始化 |  |
| `ids12` | `cfp_numerators_col_12` | DATA 初始化 |  |
| `idv12` | `cfp_denominators_col_12` | DATA 初始化 |  |
| `ids13` | `cfp_numerators_col_13` | DATA 初始化 |  |
| `idv13` | `cfp_denominators_col_13` | DATA 初始化 |  |
| `ids14` | `cfp_numerators_col_14` | DATA 初始化 |  |
| `ids15` | `cfp_numerators_col_15` | DATA 初始化 |  |
| `idv15` | `cfp_denominators_col_15` | DATA 初始化 |  |
| `ids16` | `cfp_numerators_col_16` | DATA 初始化 |  |
| `idv16` | `cfp_denominators_col_16` | DATA 初始化 |  |
| `ids17` | `cfp_numerators_col_17` | DATA 初始化 |  |
| `idv17` | `cfp_denominators_col_17` | DATA 初始化 |  |
| `ids18` | `cfp_numerators_col_18` | DATA 初始化 |  |
| `idv18` | `cfp_denominators_col_18` | DATA 初始化 |  |

### [rmew1jj](../../src/lib/spin_angular/jj_reduced_w_j_one_half.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `coef` | `reduced_w_element` | 形参 |  |
| `i10` | `w_numerators_1_0` | DATA 初始化 |  |
| `i01` | `w_numerators_0_1` | DATA 初始化 |  |

### [rmew3jj](../../src/lib/spin_angular/jj_reduced_w_j_three_halves.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `coef` | `reduced_w_element` | 形参 |  |
| `ji1` | `bra_table_index` | 局部量 |  |
| `i00` | `w_numerators_0_0` | DATA 初始化 |  |
| `i10` | `w_numerators_1_0` | DATA 初始化 |  |
| `i01` | `w_numerators_0_1` | DATA 初始化 |  |

### [rmew5jj](../../src/lib/spin_angular/jj_reduced_w_j_five_halves.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `coef` | `reduced_w_element` | 形参 |  |
| `ji1` | `bra_table_index` | 局部量 |  |
| `ji2` | `ket_table_index` | 局部量 |  |
| `i00` | `w_numerators_0_0` | DATA 初始化 |  |
| `i01` | `w_numerators_0_1` | DATA 初始化 |  |
| `i12p` | `w_numerators_1_2_group_1` | DATA 初始化 |  |
| `i12a` | `w_numerators_1_2_group_2` | DATA 初始化 |  |
| `i03p` | `w_numerators_0_3_group_1` | DATA 初始化 |  |
| `i03a` | `w_numerators_0_3_group_2` | DATA 初始化 |  |
| `i14p` | `w_numerators_1_4_group_1` | DATA 初始化 |  |
| `i14a` | `w_numerators_1_4_group_2` | DATA 初始化 |  |
| `i05p` | `w_numerators_0_5_group_1` | DATA 初始化 |  |
| `i05a` | `w_numerators_0_5_group_2` | DATA 初始化 |  |

### [rmew7bjj](../../src/lib/spin_angular/jj_reduced_w_j_seven_halves_higher_ranks.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `coef` | `reduced_w_element` | 形参 |  |
| `ll` | `two_orbital_j` | 局部量 |  |
| `iq1` | `two_bra_q` | 局部量 |  |
| `lv1` | `two_bra_seniority` | 局部量 |  |
| `il1` | `two_bra_term_j` | 局部量 |  |
| `iq2` | `two_ket_q` | 局部量 |  |
| `lv2` | `two_ket_seniority` | 局部量 |  |
| `il2` | `two_ket_term_j` | 局部量 |  |
| `ji1` | `lower_table_index` | 局部量 |  |
| `ji2` | `upper_table_index` | 局部量 |  |
| `ifaz` | `twice_phase_exponent` | 局部量 |  |
| `j` | `packed_pair_index` | 局部量 |  |
| `l` | `term_group` | 局部量 |  |
| `ipr1` | `first_group_row_offsets` | DATA 初始化 |  |
| `ipr2` | `second_group_row_offsets` | DATA 初始化 |  |
| `i12ps` | `w_numerators_1_2_group_1` | DATA 初始化 |  |
| `i12pv` | `w_denominators_1_2_group_1` | DATA 初始化 |  |
| `i03ps` | `w_numerators_0_3_group_1` | DATA 初始化 |  |
| `i03pv` | `w_denominators_0_3_group_1` | DATA 初始化 |  |
| `i14ps` | `w_numerators_1_4_group_1` | DATA 初始化 |  |
| `i14pv` | `w_denominators_1_4_group_1` | DATA 初始化 |  |
| `i05ps` | `w_numerators_0_5_group_1` | DATA 初始化 |  |
| `i05pv` | `w_denominators_0_5_group_1` | DATA 初始化 |  |
| `i16ps` | `w_numerators_1_6_group_1` | DATA 初始化 |  |
| `i16pv` | `w_denominators_1_6_group_1` | DATA 初始化 |  |
| `i07ps` | `w_numerators_0_7_group_1` | DATA 初始化 |  |
| `i07pv` | `w_denominators_0_7_group_1` | DATA 初始化 |  |
| `i12as` | `w_numerators_1_2_group_2` | DATA 初始化 |  |
| `i12av` | `w_denominators_1_2_group_2` | DATA 初始化 |  |
| `i03as` | `w_numerators_0_3_group_2` | DATA 初始化 |  |
| `i03av` | `w_denominators_0_3_group_2` | DATA 初始化 |  |
| `i14as` | `w_numerators_1_4_group_2` | DATA 初始化 |  |
| `i14av` | `w_denominators_1_4_group_2` | DATA 初始化 |  |
| `i05as` | `w_numerators_0_5_group_2` | DATA 初始化 |  |
| `i05av` | `w_denominators_0_5_group_2` | DATA 初始化 |  |
| `i16as` | `w_numerators_1_6_group_2` | DATA 初始化 |  |
| `i16av` | `w_denominators_1_6_group_2` | DATA 初始化 |  |
| `i07as` | `w_numerators_0_7_group_2` | DATA 初始化 |  |
| `i07av` | `w_denominators_0_7_group_2` | DATA 初始化 |  |

### [rmew7jj](../../src/lib/spin_angular/jj_reduced_w_j_seven_halves.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `coef` | `reduced_w_element` | 形参 |  |
| `i00` | `w_numerators_0_0` | DATA 初始化 |  |
| `i10` | `w_numerators_1_0` | DATA 初始化 |  |
| `i01` | `w_numerators_0_1` | DATA 初始化 |  |

### [rumtjj](../../src/lib/spin_angular/jj_subshell_term_decoding.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jj` | `two_orbital_j` | 形参 |  |
| `knt` | `term_index` | 形参 |  |
| `lq` | `two_q` | 形参 |  |
| `lv` | `two_seniority` | 形参 |  |
| `l` | `two_term_j` | 形参 |  |
| `kt` | `packed_term_code` | 局部量 |  |
| `kntmin` | `unoffset_term_index` | 局部量 |  |

### [rwjj](../../src/lib/spin_angular/jj_reduced_w.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j` | `two_orbital_j` | 形参 |  |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `coef` | `reduced_w_element` | 形参 |  |

### [sixj](../../src/lib/spin_angular/wigner_6j.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `two_j1` | `two_j_1` | 形参 | Wigner 6j 的第 1 个二倍角动量参数。 |
| `two_j2` | `two_j_2` | 形参 | Wigner 6j 的第 2 个二倍角动量参数。 |
| `two_j3` | `two_j_3` | 形参 | Wigner 6j 的第 3 个二倍角动量参数。 |
| `two_j4` | `two_j_4` | 形参 | Wigner 6j 的第 4 个二倍角动量参数。 |
| `two_j5` | `two_j_5` | 形参 | Wigner 6j 的第 5 个二倍角动量参数。 |
| `two_j6` | `two_j_6` | 形参 | Wigner 6j 的第 6 个二倍角动量参数。 |
| `check_triangle_conditions` | `check_triangle_conditions` | 形参 | 沿用现名。 |
| `six_j_value` | `six_j_value` | 形参 | 沿用现名。 |
| `ifa` | `twice_phase_exponent` | 局部量 |  |
| `small_six_j_cache` | `small_six_j_cache` | DATA 初始化 | 以六个 2j 值为下标的可写缓存，有 SAVE 及 OpenMP THREADPRIVATE。 |
| `uncached_six_j_sentinel` | `uncached_six_j_sentinel` | DATA 初始化 | DATA 初始化的缓存哨兵，数值和存储属性沿用原定义。 |
| `a` | `normalization_squared` | 局部量 |  |
| `store_in_six_j_cache` | `store_in_six_j_cache` | 局部量 | 沿用现名。 |

### [sixj1](../../src/lib/spin_angular/wigner_6j_with_j_one.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `i` | `two_j_1` | 形参 |  |
| `j` | `two_j_2` | 形参 |  |
| `k` | `two_j_3` | 形参 |  |
| `l` | `two_j_4` | 形参 |  |
| `m` | `two_j_5` | 形参 |  |
| `itik` | `check_triangle_conditions` | 形参 |  |
| `si` | `six_j_value` | 形参 |  |
| `ifa` | `triangle_sum` | 局部量 |  |
| `as` | `real_triangle_sum` | 局部量 |  |
| `aka` | `phase_sign` | 局部量 |  |
| `a` | `two_j_3_real` | 局部量 |  |
| `b` | `two_j_2_real` | 局部量 |  |
| `c` | `two_j_1_real` | 局部量 |  |

### [sixj2](../../src/lib/spin_angular/wigner_6j_with_j_two.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j` | `two_j_1` | 形参 |  |
| `k` | `two_j_2` | 形参 |  |
| `l` | `two_j_3` | 形参 |  |
| `m` | `two_j_4` | 形参 |  |
| `n` | `two_j_5` | 形参 |  |
| `itik` | `check_triangle_conditions` | 形参 |  |
| `si` | `six_j_value` | 形参 |  |
| `i1` | `triangle_sum` | 局部量 |  |
| `as` | `real_triangle_sum` | 局部量 |  |
| `a` | `two_j_3_real` | 局部量 |  |
| `b` | `two_j_1_real` | 局部量 |  |
| `c` | `two_j_2_real` | 局部量 |  |
| `aka` | `phase_sign` | 局部量 |  |
| `x1` | `first_numerator_term` | 局部量 |  |
| `x2` | `second_numerator_term` | 局部量 |  |
| `x3` | `third_numerator_term` | 局部量 |  |

### [sixj3](../../src/lib/spin_angular/wigner_6j_with_j_three.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j` | `two_j_1` | 形参 |  |
| `k` | `two_j_2` | 形参 |  |
| `l` | `two_j_3` | 形参 |  |
| `m` | `two_j_4` | 形参 |  |
| `n` | `two_j_5` | 形参 |  |
| `itik` | `check_triangle_conditions` | 形参 |  |
| `si` | `six_j_value` | 形参 |  |
| `i1` | `triangle_sum` | 局部量 |  |
| `as` | `real_triangle_sum` | 局部量 |  |
| `a` | `two_j_3_real` | 局部量 |  |
| `b` | `two_j_1_real` | 局部量 |  |
| `c` | `two_j_2_real` | 局部量 |  |
| `aka` | `phase_sign` | 局部量 |  |

### [sixj35](../../src/lib/spin_angular/wigner_6j_with_j_three_halves.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j` | `two_j_1` | 形参 |  |
| `k` | `two_j_2` | 形参 |  |
| `l` | `two_j_3` | 形参 |  |
| `m` | `two_j_4` | 形参 |  |
| `n` | `two_j_5` | 形参 |  |
| `itik` | `check_triangle_conditions` | 形参 |  |
| `si` | `six_j_value` | 形参 |  |
| `i1` | `triangle_sum` | 局部量 |  |
| `as` | `real_triangle_sum` | 局部量 |  |
| `a` | `two_j_3_real` | 局部量 |  |
| `b` | `two_j_1_real` | 局部量 |  |
| `c` | `two_j_2_real` | 局部量 |  |
| `aka` | `phase_sign` | 局部量 |  |

### [sixj4](../../src/lib/spin_angular/wigner_6j_with_j_four.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `jc` | `two_j_1` | 形参 |  |
| `je` | `two_j_2` | 形参 |  |
| `jd` | `two_j_3` | 形参 |  |
| `jb` | `two_j_4` | 形参 |  |
| `jf` | `two_j_5` | 形参 |  |
| `itik` | `check_triangle_conditions` | 形参 |  |
| `si` | `six_j_value` | 形参 |  |
| `a` | `j6_for_rank_3` | 局部量 |  |
| `c` | `j1` | 局部量 |  |
| `e` | `j2` | 局部量 |  |
| `d` | `j3` | 局部量 |  |
| `b` | `j4` | 局部量 |  |
| `f` | `j5` | 局部量 |  |
| `x1` | `rank_4_weight` | 局部量 |  |
| `x2` | `rank_2_weight` | 局部量 |  |
| `x3` | `rank_3_weight` | 局部量 |  |
| `s2` | `rank_2_six_j` | 局部量 |  |
| `s3` | `rank_3_six_j` | 局部量 |  |

### [sixj5](../../src/lib/spin_angular/wigner_6j_with_j_one_half.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `j` | `two_j_1` | 形参 |  |
| `k` | `two_j_2` | 形参 |  |
| `l` | `two_j_3` | 形参 |  |
| `m` | `two_j_4` | 形参 |  |
| `n` | `two_j_5` | 形参 |  |
| `itik` | `check_triangle_conditions` | 形参 |  |
| `si` | `six_j_value` | 形参 |  |
| `i1` | `triangle_sum` | 局部量 |  |
| `as` | `real_triangle_sum` | 局部量 |  |
| `a` | `two_j_3_real` | 局部量 |  |
| `b` | `two_j_2_real` | 局部量 |  |
| `c` | `two_j_1_real` | 局部量 |  |
| `aka` | `phase_sign` | 局部量 |  |

### [snrc](../../src/lib/spin_angular/breit_radial_rank_ranges.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `nd1` | `first_direct_nu` | 形参 |  |
| `nd2` | `direct_nu_count` | 形参 | 最终为直接项允许的 nu 数量。计算中暂存第一对轨道的 nu 上界。 |
| `ne1` | `first_exchange_nu` | 形参 |  |
| `ne2` | `exchange_nu_count` | 形参 | 最终为交换项允许的 nu 数量。计算中暂存第一对轨道的 nu 上界。 |
| `ibrd` | `direct_integral_type` | 形参 |  |
| `ibre` | `exchange_integral_type` | 形参 |  |
| `is` | `integral_orbitals` | 形参 |  |
| `kaps` | `two_kappas` | 形参 |  |
| `ks` | `two_j_plus_one` | 形参 |  |
| `iac` | `relative_kappa_sign_1` | 局部量 |  |
| `iad` | `relative_kappa_sign_2` | 局部量 |  |
| `nd1a` | `second_direct_lower_nu` | 局部量 |  |
| `nd2a` | `second_direct_upper_nu` | 局部量 |  |
| `ne1a` | `second_exchange_lower_nu` | 局部量 |  |
| `ne2a` | `second_exchange_upper_nu` | 局部量 |  |

### [suwjj](../../src/lib/spin_angular/jj_reduced_w_from_cfp.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `quasispin_rank` | 形参 |  |
| `k2` | `angular_rank` | 形参 |  |
| `j1` | `bra_term_index` | 形参 |  |
| `j2` | `ket_term_index` | 形参 |  |
| `ll` | `two_orbital_j` | 形参 |  |
| `suw` | `reduced_w_element` | 形参 |  |
| `kk1` | `two_quasispin_rank` | 局部量 |  |
| `kk2` | `two_angular_rank` | 局部量 |  |
| `i` | `intermediate_term_index` | 局部量 |  |
| `ip` | `first_intermediate_term` | 局部量 |  |
| `ig` | `last_intermediate_term` | 局部量 |  |
| `l2q1` | `two_bra_q` | 局部量 |  |
| `l2v1` | `two_bra_seniority` | 局部量 |  |
| `l2j1` | `two_bra_term_j` | 局部量 |  |
| `l2q2` | `two_ket_q` | 局部量 |  |
| `l2v2` | `two_ket_seniority` | 局部量 |  |
| `l2j2` | `two_ket_term_j` | 局部量 |  |
| `l2qi` | `two_intermediate_q` | 局部量 |  |
| `l2vi` | `two_intermediate_seniority` | 局部量 |  |
| `l2ji` | `two_intermediate_term_j` | 局部量 |  |
| `coef1` | `bra_reduced_cfp` | 局部量 |  |
| `coef2` | `ket_reduced_cfp` | 局部量 |  |
| `s` | `sum_over_terms` | 局部量 |  |
| `si1` | `angular_six_j` | 局部量 |  |
| `si2` | `quasispin_six_j` | 局部量 |  |

### [w1jjg](../../src/lib/spin_angular/jj_w_matrix_element_high_j.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `ww` | `w_matrix_element` | 形参 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `iq` | `occupation_change_2` | 局部量 |  |
| `iqm` | `two_abs_intermediate_projection` | 局部量 |  |
| `ie` | `twice_phase_exponent` | 局部量 |  |
| `it` | `intermediate_term_index` | 局部量 |  |
| `itk` | `bra_term_lookup_index` | 局部量 |  |
| `itd` | `ket_term_lookup_index` | 局部量 |  |
| `itp` | `first_intermediate_term` | 局部量 |  |
| `itp1` | `first_bra_range_term` | 局部量 |  |
| `itg` | `last_intermediate_term` | 局部量 |  |
| `itg1` | `last_bra_range_term` | 局部量 |  |
| `ibtt` | `two_intermediate_seniority` | 局部量 |  |
| `ibt` | `intermediate_shell_codes` | 局部量 |  |
| `enqp` | `sum_over_terms` | 局部量 |  |
| `d1` | `bra_a_matrix_element` | 局部量 |  |
| `si1` | `six_j_value` | 局部量 |  |
| `w` | `ket_a_matrix_element` | 局部量 |  |
| `bt` | `intermediate_shell_values` | 局部量 |  |

### [wap1](../../src/lib/spin_angular/jj_wa_matrix_element.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `bk2` | `total_angular_rank` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `wa` | `wa_matrix_element` | 形参 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `kk2` | `two_total_rank` | 局部量 |  |
| `ie` | `twice_phase_exponent` | 局部量 |  |
| `iq3` | `occupation_change_3` | 局部量 |  |
| `iqm` | `two_abs_intermediate_projection` | 局部量 |  |
| `it` | `intermediate_term_index` | 局部量 |  |
| `itp` | `first_intermediate_term` | 局部量 |  |
| `itg` | `last_intermediate_term` | 局部量 |  |
| `ibtt` | `two_intermediate_seniority` | 局部量 |  |
| `ibt` | `intermediate_shell_codes` | 局部量 |  |
| `enqp` | `sum_over_terms` | 局部量 |  |
| `d1` | `a_matrix_element_term` | 局部量 |  |
| `s` | `reduced_cfp` | 局部量 |  |
| `si` | `six_j_value` | 局部量 |  |
| `w` | `w_matrix_element` | 局部量 |  |
| `bt` | `intermediate_shell_values` | 局部量 |  |

### [wap1jjg](../../src/lib/spin_angular/jj_wa_matrix_element_high_j.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k1` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `bk2` | `total_angular_rank` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `wa` | `wa_matrix_element` | 形参 |  |
| `kk1` | `two_w_rank` | 局部量 |  |
| `kk2` | `two_total_rank` | 局部量 |  |
| `ie` | `twice_phase_exponent` | 局部量 |  |
| `iq3` | `occupation_change_3` | 局部量 |  |
| `iqm` | `two_abs_intermediate_projection` | 局部量 |  |
| `it` | `intermediate_term_index` | 局部量 |  |
| `itp` | `first_intermediate_term` | 局部量 |  |
| `itg` | `last_intermediate_term` | 局部量 |  |
| `ibtt` | `two_intermediate_seniority` | 局部量 |  |
| `ibt` | `intermediate_shell_codes` | 局部量 |  |
| `enqp` | `sum_over_terms` | 局部量 |  |
| `d1` | `a_matrix_element_term` | 局部量 |  |
| `si1` | `six_j_value` | 局部量 |  |
| `w` | `w_matrix_element` | 局部量 |  |
| `bt` | `intermediate_shell_values` | 局部量 |  |

### [wj1](../../src/lib/spin_angular/jj_w_matrix_element.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k2` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `wj` | `w_matrix_element` | 形参 |  |
| `k1` | `w_quasispin_rank` | 局部量 |  |
| `iq` | `net_q_projection` | 局部量 |  |
| `iqm2` | `occupation_change_2` | 局部量 |  |
| `a` | `matrix_element_factor` | 局部量 |  |
| `qq` | `combined_q_projection` | 局部量 |  |
| `w` | `reduced_w_element` | 局部量 |  |
| `wk1` | `real_quasispin_rank` | 局部量 |  |

### [ww1](../../src/lib/spin_angular/jj_ww_scalar_matrix_element.f90)

| 原名 | 迁移后 | 角色 | 备注 |
| --- | --- | --- | --- |
| `k2` | `w_angular_rank` | 形参 |  |
| `ik` | `bra_shell_codes` | 形参 |  |
| `id` | `ket_shell_codes` | 形参 |  |
| `qm1` | `q_projection_1` | 形参 |  |
| `qm2` | `q_projection_2` | 形参 |  |
| `qm3` | `q_projection_3` | 形参 |  |
| `qm4` | `q_projection_4` | 形参 |  |
| `bk` | `bra_shell_values` | 形参 |  |
| `bd` | `ket_shell_values` | 形参 |  |
| `ww` | `ww_scalar_matrix_element` | 形参 |  |
| `kk2` | `two_w_rank` | 局部量 |  |
| `kk6` | `two_total_rank` | 局部量 |  |
| `iq` | `intermediate_occupation_change` | 局部量 |  |
| `iqm` | `two_abs_intermediate_projection` | 局部量 |  |
| `iq3` | `occupation_change_3` | 局部量 |  |
| `iq4` | `occupation_change_4` | 局部量 |  |
| `ie` | `twice_phase_exponent` | 局部量 |  |
| `ie1` | `partial_phase_sum` | 局部量 |  |
| `it` | `intermediate_term_index` | 局部量 |  |
| `itp` | `first_intermediate_term` | 局部量 |  |
| `itg` | `last_intermediate_term` | 局部量 |  |
| `ibtt` | `two_intermediate_seniority` | 局部量 |  |
| `ibt` | `intermediate_shell_codes` | 局部量 |  |
| `enqp` | `sum_over_terms` | 局部量 |  |
| `d1` | `bra_w_matrix_element` | 局部量 |  |
| `w` | `ket_w_matrix_element` | 局部量 |  |
| `bt` | `intermediate_shell_values` | 局部量 |  |
