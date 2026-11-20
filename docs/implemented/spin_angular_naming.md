# 自旋–角向库命名迁移

日期 2026-11-20。基准为 suite 的 `2d4c721` 及同一工作树中已实施的变量与共享库迁移。

| 项目 | 旧名 | 现名 |
| --- | --- | --- |
| 源码目录 | `src/lib/librang90` | [`src/lib/spin_angular`](../../src/lib/spin_angular/CMakeLists.txt) |
| CMake 目标 | `rang90` | `spin_angular` |
| 静态库 | `librang90.a` | `libspin_angular.a` |

89 个源文件采用小写下划线命名，文件基名就是其中定义的 Fortran 模块名。
名称按职责区分 Wigner 与 Racah 系数、准自旋及子壳层项、张量矩阵元、重耦合和单／双电子角向系数。
`A` 与 `W` 沿用角向代数公式中的张量记号。`four_shell_transfer_13_24` 等名称的数字指相互配对的子壳层编号。

仓库内的 `use` 语句、库链接、F2018 必需目标清单和 `jj2lsj` 的 OpenMP 源文件清单采用新名称。
过程名和实参接口沿用原有定义。
外部 Fortran 调用方需要更新 `use` 和库链接名称，并重新编译。
形参和局部量见[变量与常量对照](librang90_variable_migration.md)。
角向计算公式见[自旋–角向系数说明](../../src/lib/spin_angular/readme.md)。

## 源文件与模块对照

| 原文件 | 原模块 | 现文件 | 现模块 | 职责 |
| --- | --- | --- | --- | --- |
| `Gracah1.f90` | `gracah1_i` | [racah_coefficients_scott_fischer.f90](../../src/lib/spin_angular/racah_coefficients_scott_fischer.f90) | `racah_coefficients_scott_fischer` | Scott–Froese Fischer 版本的 Racah 系数。 |
| `dracah.f90` | `dracah_i` | [racah_coefficients.f90](../../src/lib/spin_angular/racah_coefficients.f90) | `racah_coefficients` | Wigner 6j 通用路径使用的 Racah 系数。 |
| `sixj.f90` | `sixj_i` | [wigner_6j.f90](../../src/lib/spin_angular/wigner_6j.f90) | `wigner_6j` | Wigner 6j 符号的选择规则、缓存与求值入口。 |
| `sixj1.f90` | `sixj1_i` | [wigner_6j_with_j_one.f90](../../src/lib/spin_angular/wigner_6j_with_j_one.f90) | `wigner_6j_with_j_one` | 含 j = 1 的解析特例。 |
| `sixj2.f90` | `sixj2_i` | [wigner_6j_with_j_two.f90](../../src/lib/spin_angular/wigner_6j_with_j_two.f90) | `wigner_6j_with_j_two` | 含 j = 2 的解析特例。 |
| `sixj3.f90` | `sixj3_i` | [wigner_6j_with_j_three.f90](../../src/lib/spin_angular/wigner_6j_with_j_three.f90) | `wigner_6j_with_j_three` | 含 j = 3 的解析特例。 |
| `sixj35.f90` | `sixj35_i` | [wigner_6j_with_j_three_halves.f90](../../src/lib/spin_angular/wigner_6j_with_j_three_halves.f90) | `wigner_6j_with_j_three_halves` | 含 j = 3/2 的解析特例。 |
| `sixj4.f90` | `sixj4_i` | [wigner_6j_with_j_four.f90](../../src/lib/spin_angular/wigner_6j_with_j_four.f90) | `wigner_6j_with_j_four` | 含 j = 4 的递推特例。 |
| `sixj5.f90` | `sixj5_i` | [wigner_6j_with_j_one_half.f90](../../src/lib/spin_angular/wigner_6j_with_j_one_half.f90) | `wigner_6j_with_j_one_half` | 含 j = 1/2 的解析特例。 |
| `nine.f90` | `nine_i` | [wigner_9j.f90](../../src/lib/spin_angular/wigner_9j.f90) | `wigner_9j` | Wigner 9j 符号的三角条件与 6j 乘积求和。 |
| `nine0.f90` | `nine0_i` | [wigner_9j_with_zero.f90](../../src/lib/spin_angular/wigner_9j_with_zero.f90) | `wigner_9j_with_zero` | 含零角动量时将 9j 化为单个 6j。 |
| `ixjtik.f90` | `ixjtik_i` | [wigner_6j_triangle_selection.f90](../../src/lib/spin_angular/wigner_6j_triangle_selection.f90) | `wigner_6j_triangle_selection` | 检查 6j 的四组三角条件。 |
| `ittk.f90` | `ittk_i` | [angular_triangle_and_parity.f90](../../src/lib/spin_angular/angular_triangle_and_parity.f90) | `angular_triangle_and_parity` | 检查三角不等式与三量之和的奇偶性。 |
| `itrexg.f90` | `itrexg_i` | [shared_tensor_rank_range.f90](../../src/lib/spin_angular/shared_tensor_rank_range.f90) | `shared_tensor_rank_range` | 求两组轨道角动量同时允许的二倍张量阶数区间。 |
| `coulom.f90` | `coulom_i` | [coulomb_angular_factor.f90](../../src/lib/spin_angular/coulomb_angular_factor.f90) | `coulomb_angular_factor` | 由两个球张量约化元组成 Coulomb 角向因子。 |
| `snrc.f90` | `snrc_i` | [breit_radial_rank_ranges.f90](../../src/lib/spin_angular/breit_radial_rank_ranges.f90) | `breit_radial_rank_ranges` | 确定 Breit 直接项、交换项的径向阶数范围与积分类型。 |
| `eile.f90` | `eile_i` | [three_subshell_ordering.f90](../../src/lib/spin_angular/three_subshell_ordering.f90) | `three_subshell_ordering` | 排序三个互异的子壳层位置。 |
| `jfaze.f90` | `jfaze_i` | [subshell_permutation_sign.f90](../../src/lib/spin_angular/subshell_permutation_sign.f90) | `subshell_permutation_sign` | 由四个子壳层位置的逆序数确定置换符号。 |
| `jthn.f90` | `jthn_i` | [packed_integer_digit.f90](../../src/lib/spin_angular/packed_integer_digit.f90) | `packed_integer_digit` | 读取指定基数表示的整数位。 |
| `mes.f90` | `mes_i` | [subshell_tensor_diagnostics.f90](../../src/lib/spin_angular/subshell_tensor_diagnostics.f90) | `subshell_tensor_diagnostics` | 报告子壳层张量的项编码或占据数错误。 |
| `izas1.f90` | `izas1_i` | [quasispin_projection_selection.f90](../../src/lib/spin_angular/quasispin_projection_selection.f90) | `quasispin_projection_selection` | 检查两侧准自旋投影的范围和奇偶。 |
| `c0t5s.f90` | `c0t5s_i` | [quasispin_cg_rank_half.f90](../../src/lib/spin_angular/quasispin_cg_rank_half.f90) | `quasispin_cg_rank_half` | 准自旋张量秩为 1/2 的 Clebsch–Gordan 系数。 |
| `c1e0sm.f90` | `c1e0sm_i` | [quasispin_cg_rank_one_zero_projection.f90](../../src/lib/spin_angular/quasispin_cg_rank_one_zero_projection.f90) | `quasispin_cg_rank_one_zero_projection` | 秩为 1、投影为 0 的准自旋 Clebsch–Gordan 系数。 |
| `c1e1sm.f90` | `c1e1sm_i` | [quasispin_cg_rank_one_unit_projection.f90](../../src/lib/spin_angular/quasispin_cg_rank_one_unit_projection.f90) | `quasispin_cg_rank_one_unit_projection` | 秩为 1、投影为 ±1 的准自旋 Clebsch–Gordan 系数。 |
| `cle0sm.f90` | `cle0sm_i` | [quasispin_cg_zero_projection.f90](../../src/lib/spin_angular/quasispin_cg_zero_projection.f90) | `quasispin_cg_zero_projection` | 秩为 0 或 1、投影为 0 的准自旋 Clebsch–Gordan 系数。 |
| `nmtejj.f90` | `nmtejj_i` | [jj_subshell_term_lookup.f90](../../src/lib/spin_angular/jj_subshell_term_lookup.f90) | `jj_subshell_term_lookup` | 从 j、Q、J 与占据数查找子壳层项索引。 |
| `rumtjj.f90` | `rumtjj_i` | [jj_subshell_term_decoding.f90](../../src/lib/spin_angular/jj_subshell_term_decoding.f90) | `jj_subshell_term_decoding` | 从项索引解码二倍 Q、seniority 和 J。 |
| `itjj.f90` | `itjj_i` | [jj_intermediate_terms_same_parity.f90](../../src/lib/spin_angular/jj_intermediate_terms_same_parity.f90) | `jj_intermediate_terms_same_parity` | 检查同奇偶占据的 bra 和 ket，给出中间项范围。 |
| `itjj2.f90` | `itjj2_i` | [jj_intermediate_terms_ket_parity.f90](../../src/lib/spin_angular/jj_intermediate_terms_ket_parity.f90) | `jj_intermediate_terms_ket_parity` | 给出与 ket 侧占据奇偶相同的中间项范围。 |
| `itjj3.f90` | `itjj3_i` | [jj_intermediate_terms_bra_parity.f90](../../src/lib/spin_angular/jj_intermediate_terms_bra_parity.f90) | `jj_intermediate_terms_bra_parity` | 给出与 bra 侧占据奇偶相同的中间项范围。 |
| `perko1.f90` | `perko1_i` | [jj_subshell_state_packing.f90](../../src/lib/spin_angular/jj_subshell_state_packing.f90) | `jj_subshell_state_packing` | 将一对 CSF 中某个子壳层的量子数写入整数和实数槽。 |
| `perko2.f90` | `perko2_i` | [jj_active_subshell_states.f90](../../src/lib/spin_angular/jj_active_subshell_states.f90) | `jj_active_subshell_states` | 填充至多四个活跃子壳层的张量工作区。 |
| `Rmeajj.f90` | `rmeajj_i` | [jj_reduced_cfp.f90](../../src/lib/spin_angular/jj_reduced_cfp.f90) | `jj_reduced_cfp` | 单个产生或湮灭张量的完全约化分数亲代系数入口。 |
| `Rmeajj9.f90` | `rmeajj9_i` | [jj_reduced_cfp_j_nine_halves.f90](../../src/lib/spin_angular/jj_reduced_cfp_j_nine_halves.f90) | `jj_reduced_cfp_j_nine_halves` | j = 9/2 子壳层的完全约化 CFP 表。 |
| `Rmeajj11.f90` | `rmeajj11_i` | [jj_reduced_cfp_two_electron_terms.f90](../../src/lib/spin_angular/jj_reduced_cfp_two_electron_terms.f90) | `jj_reduced_cfp_two_electron_terms` | 从双电子项与单电子项计算高 j 子壳层的完全约化 CFP。 |
| `Rwjj.f90` | `rwjj_i` | [jj_reduced_w.f90](../../src/lib/spin_angular/jj_reduced_w.f90) | `jj_reduced_w` | 按单电子 j 选择双张量 W 的完全约化矩阵元。 |
| `Rmew1jj.f90` | `rmew1jj_i` | [jj_reduced_w_j_one_half.f90](../../src/lib/spin_angular/jj_reduced_w_j_one_half.f90) | `jj_reduced_w_j_one_half` | j = 1/2 的完全约化 W 表。 |
| `Rmew3jj.f90` | `rmew3jj_i` | [jj_reduced_w_j_three_halves.f90](../../src/lib/spin_angular/jj_reduced_w_j_three_halves.f90) | `jj_reduced_w_j_three_halves` | j = 3/2 的完全约化 W 表。 |
| `Rmew5jj.f90` | `rmew5jj_i` | [jj_reduced_w_j_five_halves.f90](../../src/lib/spin_angular/jj_reduced_w_j_five_halves.f90) | `jj_reduced_w_j_five_halves` | j = 5/2 的完全约化 W 表。 |
| `Rmew7jj.f90` | `rmew7jj_i` | [jj_reduced_w_j_seven_halves.f90](../../src/lib/spin_angular/jj_reduced_w_j_seven_halves.f90) | `jj_reduced_w_j_seven_halves` | j = 7/2 的完全约化 W 表入口及低秩条目。 |
| `Rmew7bjj.f90` | `rmew7bjj_i` | [jj_reduced_w_j_seven_halves_higher_ranks.f90](../../src/lib/spin_angular/jj_reduced_w_j_seven_halves_higher_ranks.f90) | `jj_reduced_w_j_seven_halves_higher_ranks` | j = 7/2 的高秩完全约化 W 表。 |
| `suwjj.f90` | `suwjj_i` | [jj_reduced_w_from_cfp.f90](../../src/lib/spin_angular/jj_reduced_w_from_cfp.f90) | `jj_reduced_w_from_cfp` | 对中间项的两个完全约化 CFP 求和构造 W。 |
| `a1jj.f90` | `a1jj_i` | [jj_a_matrix_element_few_electrons.f90](../../src/lib/spin_angular/jj_a_matrix_element_few_electrons.f90) | `jj_a_matrix_element_few_electrons` | 少电子高 j 路径中的单个 A 张量矩阵元。 |
| `wj1.f90` | `wj1_i` | [jj_w_matrix_element.f90](../../src/lib/spin_angular/jj_w_matrix_element.f90) | `jj_w_matrix_element` | 带准自旋投影的 W 张量矩阵元。 |
| `w1jjg.f90` | `w1jjg_i` | [jj_w_matrix_element_high_j.f90](../../src/lib/spin_angular/jj_w_matrix_element_high_j.f90) | `jj_w_matrix_element_high_j` | 高 j 的 W 矩阵元，由两个 A 张量求和。 |
| `awp1.f90` | `awp1_i` | [jj_aw_matrix_element.f90](../../src/lib/spin_angular/jj_aw_matrix_element.f90) | `jj_aw_matrix_element` | 按 A、W 顺序耦合的同壳三张量矩阵元。 |
| `awp1jjg.f90` | `awp1jjg_i` | [jj_aw_matrix_element_high_j.f90](../../src/lib/spin_angular/jj_aw_matrix_element_high_j.f90) | `jj_aw_matrix_element_high_j` | 高 j 的 AW 矩阵元。 |
| `wap1.f90` | `wap1_i` | [jj_wa_matrix_element.f90](../../src/lib/spin_angular/jj_wa_matrix_element.f90) | `jj_wa_matrix_element` | 按 W、A 顺序耦合的同壳三张量矩阵元。 |
| `wap1jjg.f90` | `wap1jjg_i` | [jj_wa_matrix_element_high_j.f90](../../src/lib/spin_angular/jj_wa_matrix_element_high_j.f90) | `jj_wa_matrix_element_high_j` | 高 j 的 WA 矩阵元。 |
| `ww1.f90` | `ww1_i` | [jj_ww_scalar_matrix_element.f90](../../src/lib/spin_angular/jj_ww_scalar_matrix_element.f90) | `jj_ww_scalar_matrix_element` | 同一子壳层的两个 W 耦合为标量。 |
| `gg12.f90` | `gg12_i` | [subshell_product_a_a.f90](../../src/lib/spin_angular/subshell_product_a_a.f90) | `subshell_product_a_a` | 两个子壳层各含一个 A 张量。 |
| `gg1122.f90` | `gg1122_i` | [subshell_product_w_w.f90](../../src/lib/spin_angular/subshell_product_w_w.f90) | `subshell_product_w_w` | 两个子壳层各含一个 W 张量。 |
| `gg1112.f90` | `gg1112_i` | [subshell_product_wa_a.f90](../../src/lib/spin_angular/subshell_product_wa_a.f90) | `subshell_product_wa_a` | 第一壳层含 WA，第二壳层含 A。 |
| `gg1222.f90` | `gg1222_i` | [subshell_product_a_aw.f90](../../src/lib/spin_angular/subshell_product_a_aw.f90) | `subshell_product_a_aw` | 第一壳层含 A，第二壳层含 AW。 |
| `gg1233.f90` | `gg1233_i` | [subshell_product_a_a_w.f90](../../src/lib/spin_angular/subshell_product_a_a_w.f90) | `subshell_product_a_a_w` | 两个子壳层各含 A，第三壳层含 W。 |
| `gg1234.f90` | `gg1234_i` | [subshell_product_a_a_a_a.f90](../../src/lib/spin_angular/subshell_product_a_a_a_a.f90) | `subshell_product_a_a_a_a` | 四个子壳层各含一个 A 张量。 |
| `diaga1.f90` | `diaga1_i` | [first_active_shell_recoupling.f90](../../src/lib/spin_angular/first_active_shell_recoupling.f90) | `first_active_shell_recoupling` | 第一个活跃子壳层处的 6j 重耦合因子。 |
| `diaga2.f90` | `diaga2_i` | [later_active_shell_recoupling.f90](../../src/lib/spin_angular/later_active_shell_recoupling.f90) | `later_active_shell_recoupling` | 后续活跃子壳层处的 6j 重耦合因子。 |
| `diaga3.f90` | `diaga3_i` | [spectator_shell_chain_recoupling.f90](../../src/lib/spin_angular/spectator_shell_chain_recoupling.f90) | `spectator_shell_chain_recoupling` | 两个活跃位置之间的旁观子壳层重耦合因子乘积。 |
| `diaga4.f90` | `diaga4_i` | [pair_tensor_recoupling.f90](../../src/lib/spin_angular/pair_tensor_recoupling.f90) | `pair_tensor_recoupling` | 两个张量耦合时的 9j 重耦合因子。 |
| `diaga5.f90` | `diaga5_i` | [final_shell_recoupling.f90](../../src/lib/spin_angular/final_shell_recoupling.f90) | `final_shell_recoupling` | 耦合链最后一个子壳层处的 6j 因子。 |
| `reco.f90` | `reco_i` | [scalar_interaction_recoupling_selection.f90](../../src/lib/spin_angular/scalar_interaction_recoupling_selection.f90) | `scalar_interaction_recoupling_selection` | 根据活跃子壳层及耦合链比较两侧 CSF，筛选标量相互作用。 |
| `reco2.f90` | `reco2_i` | [two_subshell_recoupling.f90](../../src/lib/spin_angular/two_subshell_recoupling.f90) | `two_subshell_recoupling` | 两个活跃子壳层的重耦合。 |
| `reco3.f90` | `reco3_i` | [three_subshell_recoupling.f90](../../src/lib/spin_angular/three_subshell_recoupling.f90) | `three_subshell_recoupling` | 三个已排序的活跃子壳层的重耦合。 |
| `reco4.f90` | `reco4_i` | [four_subshell_recoupling.f90](../../src/lib/spin_angular/four_subshell_recoupling.f90) | `four_subshell_recoupling` | 四个活跃子壳层的重耦合。 |
| `rec3.f90` | `rec3_i` | [permuted_three_subshell_recoupling.f90](../../src/lib/spin_angular/permuted_three_subshell_recoupling.f90) | `permuted_three_subshell_recoupling` | 排序三个活跃位置，再计入张量置换相位。 |
| `recoonescalar.f90` | `recoonescalar_i` | [one_electron_scalar_recoupling_selection.f90](../../src/lib/spin_angular/one_electron_scalar_recoupling_selection.f90) | `one_electron_scalar_recoupling_selection` | 单电子标量路径的重耦合筛选，包含有效子壳层数。 |
| `recop00.f90` | `recop00_i` | [one_electron_tensor_recoupling_selection.f90](../../src/lib/spin_angular/one_electron_tensor_recoupling_selection.f90) | `one_electron_tensor_recoupling_selection` | 一般单电子张量的三角条件与旁观子壳层筛选。 |
| `recop1.f90` | `recop1_i` | [one_electron_same_subshell_recoupling.f90](../../src/lib/spin_angular/one_electron_same_subshell_recoupling.f90) | `one_electron_same_subshell_recoupling` | 单电子算符作用于同一子壳层时的重耦合。 |
| `recop2.f90` | `recop2_i` | [one_electron_two_subshell_recoupling.f90](../../src/lib/spin_angular/one_electron_two_subshell_recoupling.f90) | `one_electron_two_subshell_recoupling` | 单电子算符连接两个子壳层时的重耦合。 |
| `onescalar.f90` | `onescalar_i` | [one_electron_scalar_coefficients.f90](../../src/lib/spin_angular/one_electron_scalar_coefficients.f90) | `one_electron_scalar_coefficients` | 单电子偶标量算符的角向系数入口。 |
| `onescalar1.f90` | `onescalar1_i` | [one_electron_scalar_same_subshell.f90](../../src/lib/spin_angular/one_electron_scalar_same_subshell.f90) | `one_electron_scalar_same_subshell` | 单电子标量算符作用于同一子壳层。 |
| `onescalar2.f90` | `onescalar2_i` | [one_electron_scalar_between_subshells.f90](../../src/lib/spin_angular/one_electron_scalar_between_subshells.f90) | `one_electron_scalar_between_subshells` | 单电子标量算符连接两个子壳层。 |
| `oneparticlejj.f90` | `oneparticlejj_i` | [one_electron_tensor_coefficients.f90](../../src/lib/spin_angular/one_electron_tensor_coefficients.f90) | `one_electron_tensor_coefficients` | 给定张量秩和宇称的单电子角向系数入口。 |
| `oneparticlejj1.f90` | `oneparticlejj1_i` | [one_electron_tensor_same_subshell.f90](../../src/lib/spin_angular/one_electron_tensor_same_subshell.f90) | `one_electron_tensor_same_subshell` | 一般单电子张量作用于同一子壳层。 |
| `oneparticlejj2.f90` | `oneparticlejj2_i` | [one_electron_tensor_between_subshells.f90](../../src/lib/spin_angular/one_electron_tensor_between_subshells.f90) | `one_electron_tensor_between_subshells` | 一般单电子张量连接两个子壳层。 |
| `rkco_gg.f90` | `rkco_gg_i` | [two_electron_angular_coefficients.f90](../../src/lib/spin_angular/two_electron_angular_coefficients.f90) | `two_electron_angular_coefficients` | 按两侧占据差分派 Coulomb 和 Breit 双电子项。 |
| `el1.f90` | `el1_i` | [two_electron_unchanged_occupations.f90](../../src/lib/spin_angular/two_electron_unchanged_occupations.f90) | `two_electron_unchanged_occupations` | 两侧各子壳层占据数相同的双电子项。 |
| `el2.f90` | `el2_i` | [two_electron_pair_transfer.f90](../../src/lib/spin_angular/two_electron_pair_transfer.f90) | `two_electron_pair_transfer` | 两个电子成对转移于两个子壳层。 |
| `el3.f90` | `el3_i` | [two_electron_single_transfer.f90](../../src/lib/spin_angular/two_electron_single_transfer.f90) | `two_electron_single_transfer` | 双电子算符中净占据变化等效于一个电子转移的分派入口。 |
| `el31.f90` | `el31_i` | [single_transfer_triple_on_first_shell.f90](../../src/lib/spin_angular/single_transfer_triple_on_first_shell.f90) | `single_transfer_triple_on_first_shell` | 净转移一个电子，四个算符中三个作用于编号为 1 的子壳层。 |
| `el32.f90` | `el32_i` | [single_transfer_triple_on_second_shell.f90](../../src/lib/spin_angular/single_transfer_triple_on_second_shell.f90) | `single_transfer_triple_on_second_shell` | 净转移一个电子，四个算符中三个作用于编号为 2 的子壳层。 |
| `el33.f90` | `el33_i` | [single_transfer_with_spectator_shell.f90](../../src/lib/spin_angular/single_transfer_with_spectator_shell.f90) | `single_transfer_with_spectator_shell` | 净转移一个电子，第三个子壳层承载两个算符且占据数相同。 |
| `el4.f90` | `el4_i` | [two_electron_three_shell_transfer.f90](../../src/lib/spin_angular/two_electron_three_shell_transfer.f90) | `two_electron_three_shell_transfer` | 两个电子在三个子壳层间转移的分派入口。 |
| `el41.f90` | `el41_i` | [three_shell_pair_split_or_join.f90](../../src/lib/spin_angular/three_shell_pair_split_or_join.f90) | `three_shell_pair_split_or_join` | 两个电子从同一壳层分至两个壳层，或由两个壳层汇入同一壳层。 |
| `el5.f90` | `el5_i` | [two_electron_four_shell_transfer.f90](../../src/lib/spin_angular/two_electron_four_shell_transfer.f90) | `two_electron_four_shell_transfer` | 四个子壳层各变化一个电子时的双电子分派入口。 |
| `el51.f90` | `el51_i` | [four_shell_transfer_12_34.f90](../../src/lib/spin_angular/four_shell_transfer_12_34.f90) | `four_shell_transfer_12_34` | 四壳层转移中，按 1–2 与 3–4 分组。 |
| `el52.f90` | `el52_i` | [four_shell_transfer_13_24.f90](../../src/lib/spin_angular/four_shell_transfer_13_24.f90) | `four_shell_transfer_13_24` | 四壳层转移中，按 1–3 与 2–4 分组。 |
| `el53.f90` | `el53_i` | [four_shell_transfer_14_23.f90](../../src/lib/spin_angular/four_shell_transfer_14_23.f90) | `four_shell_transfer_14_23` | 四壳层转移中，按 1–4 与 2–3 分组。 |

## 验证

- GNU Fortran 15.3 的 Release 构建覆盖串行、MPI 和 `jj2lsj` OpenMP。隔离安装验证了
  `libspin_angular.a` 与 89 个模块文件，独立 Fortran 程序从安装目录链接并调用了 `sixj`。
- CTest 首轮通过 107 项。工作树中的共享库迁移补齐来源清单及 QED 测试的模块引用后，
  其余 3 项复验通过。
- 同配置对照改名前后的二进制。串行、2 rank 和 4 rank 的 `rangular` 共 42 个 MCP 文件逐字节一致，
  2 rank 的 `rci_mpi` 在 Coulomb 和 Breit 两种设置下生成的混合文件也分别逐字节一致。
