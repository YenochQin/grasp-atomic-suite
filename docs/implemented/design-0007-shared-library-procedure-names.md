# GRASP 共享库过程命名

2026-11-20 在共享库的模块命名迁移基础上更新了过程名。
表中的模块路径指向过程的实现文件。外部 Fortran 调用方可按表更新 `USE, ONLY` 导入，
使用重新安装的 `.mod` 文件和静态库重新编译、链接。使用 CSF 字节写入入口时，关键字为
`unpacked_value`、`byte_index` 和 `packed_values`。

`radial_node_count` 和 `csf_byte_packing` 直接定义 `count_radial_nodes` 和 `pack_csf_subshell_byte`。
Davidson 求解入口为 `solve_davidson_eigenpairs`，矩阵–向量回调接口为 `davidson_matrix_vector_operator`。
双电子角系数入口为 `calculate_two_electron_angular_coefficients`。

2026-11-21 的语义命名检查更新了以下共享接口。外部 Fortran 调用方按新名更新 `USE`
和过程引用后重新编译。Rust 描述符读库的 C 头文件补齐了形参名称，沿用原 C ABI。

| 原过程或变量名 | 当前名称 |
| --- | --- |
| `mpi_runtime::mpix_bytes` | `create_contiguous_mpi_byte_datatype` |
| `mpi_utilities::mpix_startup` | `start_mpi_application` |
| `mpi_utilities::mpix_shutdown` | `finish_mpi_application` |
| `davidson_eigensolver::dvd_operator` | `davidson_matrix_vector_operator` |
| `debug_output_settings::cutoff` | `matrix_element_screening_cutoff` |

## 核心库

| 原过程名 | 当前过程名 | 源文件 |
| --- | --- | --- |
| `alcbuf` | `manage_angular_coefficient_buffer` | [angular_coefficient_buffer](../../src/lib/core/angular_coefficient_buffer.f90) |
| `arctan` | `nonnegative_angle` | [nonnegative_atan2](../../src/lib/core/nonnegative_atan2.f90) |
| `calen` | `read_calendar_clock` | [calendar_clock](../../src/lib/core/calendar_clock.f90) |
| `cgamma` | `complex_gamma_value` | [complex_gamma](../../src/lib/core/complex_gamma.f90) |
| `clrx` | `relativistic_three_j_value` | [relativistic_three_j](../../src/lib/core/relativistic_three_j.f90) |
| `convrt` | `format_integer` | [integer_format](../../src/lib/core/integer_format.f90) |
| `convrt2` | `format_integer_with_context` | [integer_format_with_context](../../src/lib/core/integer_format_with_context.f90) |
| `convrt_double` | `format_half_integer` | [half_integer_format](../../src/lib/core/half_integer_format.f90) |
| `cord` | `calculate_closed_shell_mcp` | [closed_shell_mcp](../../src/lib/core/closed_shell_mcp.f90) |
| `count`（实现 `count_impl`） | `count_radial_nodes` | [radial_node_count](../../src/lib/core/radial_node_count.f90) |
| `cre` | `reduced_angular_element` | [reduced_angular_matrix_element](../../src/lib/core/reduced_angular_matrix_element.f90) |
| `dcbsrw` | `evaluate_dirac_coulomb_orbital` | [dirac_coulomb_radial_orbital](../../src/lib/core/dirac_coulomb_radial_orbital.f90) |
| `dinit` | `initialize_real_vector` | [real_vector_initialization](../../src/lib/core/real_vector_initialization.f90) |
| `dmerge` | `sparse_dot_scatter` | [sparse_dot_and_scatter](../../src/lib/core/sparse_dot_and_scatter.f90) |
| `dpbdt` | `differentiate_radial_wavefunctions` | [radial_wavefunction_derivatives](../../src/lib/core/radial_wavefunction_derivatives.f90) |
| `draw` | `plot_radial_orbital_ascii` | [radial_ascii_plot](../../src/lib/core/radial_ascii_plot.f90) |
| `es` | `evaluate_alternating_exponential_series` | [alternating_exponential_series](../../src/lib/core/alternating_exponential_series.f90) |
| `factt` | `initialize_log_factorial_table` | [log_factorial_table](../../src/lib/core/log_factorial_table.f90) |
| `getrsl` | `select_subshell_labels` | [subshell_label_selection](../../src/lib/core/subshell_label_selection.f90) |
| `getyn` | `read_yes_no` | [yes_no_input](../../src/lib/core/yes_no_input.f90) |
| `ichkq1` | `check_one_electron_csf_pair` | [one_electron_csf_filter](../../src/lib/core/one_electron_csf_filter.f90) |
| `ichkq2` | `check_two_electron_csf_pair` | [two_electron_csf_filter](../../src/lib/core/two_electron_csf_filter.f90) |
| `ichop` | `classify_subshell_occupation` | [subshell_occupation_class](../../src/lib/core/subshell_occupation_class.f90) |
| `icopy` | `copy_integer_vector` | [integer_vector_copy](../../src/lib/core/integer_vector_copy.f90) |
| `iniest2` | `estimate_sparse_initial_basis` | [sparse_eigensolver_initial_basis](../../src/lib/core/sparse_eigensolver_initial_basis.f90) |
| `interp` | `interpolate_aitken` | [aitken_interpolation](../../src/lib/core/aitken_interpolation.f90) |
| `intrpq` | `interpolate_radial_orbital` | [radial_orbital_interpolation](../../src/lib/core/radial_orbital_interpolation.f90) |
| `iq` | `csf_occupation` | [csf_quantum_queries](../../src/lib/core/csf_quantum_queries.f90) |
| `ispar` | `csf_parity_code` | [csf_quantum_queries](../../src/lib/core/csf_quantum_queries.f90) |
| `items` | `parse_level_selection` | [level_selection_parser](../../src/lib/core/level_selection_parser.f90) |
| `itjpo` | `csf_total_j_code` | [csf_quantum_queries](../../src/lib/core/csf_quantum_queries.f90) |
| `itrig` | `satisfies_angular_triangle` | [angular_momentum_triangle](../../src/lib/core/angular_momentum_triangle.f90) |
| `jcup` | `csf_intermediate_coupling` | [csf_quantum_queries](../../src/lib/core/csf_quantum_queries.f90) |
| `jqs` | `csf_subshell_quantum_code` | [csf_quantum_queries](../../src/lib/core/csf_quantum_queries.f90) |
| `ldigit` | `is_decimal_digit` | [decimal_digit_check](../../src/lib/core/decimal_digit_check.f90) |
| `lodiso` | `read_isotope_records` | [isotope_data_records](../../src/lib/core/isotope_data_records.f90) |
| `lodrwf` | `read_radial_wavefunction_records` | [radial_wavefunction_records](../../src/lib/core/radial_wavefunction_records.f90) |
| `lodstate` | `read_block_level_selection` | [block_level_selection](../../src/lib/core/block_level_selection.f90) |
| `ltab` | `lookup_subshell_parents` | [subshell_parent_table](../../src/lib/core/subshell_parent_table.f90) |
| `nucpot` | `calculate_nuclear_coulomb_potential` | [nuclear_coulomb_potential](../../src/lib/core/nuclear_coulomb_potential.f90) |
| `openfl` | `open_calculation_file` | [file_open](../../src/lib/core/file_open.f90) |
| `orthsc` | `orthogonalize_radial_orbitals` | [radial_orbital_orthogonalization](../../src/lib/core/radial_orbital_orthogonalization.f90) |
| `pack`（实现 `pack_impl`） | `pack_csf_subshell_byte` | [csf_byte_packing](../../src/lib/core/csf_byte_packing.f90) |
| `posfile` | `position_calculation_file` | [file_positioning](../../src/lib/core/file_positioning.f90) |
| `quad` | `integrate_radial_grid` | [radial_quadrature](../../src/lib/core/radial_quadrature.f90) |
| `radgrd` | `initialize_radial_grid` | [radial_grid_setup](../../src/lib/core/radial_grid_setup.f90) |
| `rint` | `integrate_weighted_radial_overlap` | [weighted_radial_overlap](../../src/lib/core/weighted_radial_overlap.f90) |
| `rinti` | `integrate_one_electron_radial` | [one_electron_radial_integral](../../src/lib/core/one_electron_radial_integral.f90) |
| `setcon` | `initialize_physical_constants` | [physical_constants](../../src/lib/core/physical_constants.f90) |
| `setcsla` | `load_angular_csf_basis` | [csf_basis_loader](../../src/lib/core/csf_basis_loader.f90) |
| `setiso` | `read_isotope_data` | [isotope_data_input](../../src/lib/core/isotope_data_input.f90) |
| `setj` | `initialize_angular_recoupling_tables` | [angular_recoupling_tables](../../src/lib/core/angular_recoupling_tables.f90) |
| `setmc` | `initialize_machine_constants` | [machine_constants](../../src/lib/core/machine_constants.f90) |
| `setpot` | `calculate_radial_dirac_coefficients` | [radial_dirac_coefficients](../../src/lib/core/radial_dirac_coefficients.f90) |
| `setqic` | `initialize_radial_stencil_weights` | [radial_stencil_weights](../../src/lib/core/radial_stencil_weights.f90) |
| `setqna` | `prepare_csf_pair_quantum_numbers` | [csf_pair_quantum_numbers](../../src/lib/core/csf_pair_quantum_numbers.f90) |
| `setrwfa` | `read_radial_wavefunction_file` | [radial_wavefunction_input](../../src/lib/core/radial_wavefunction_input.f90) |
| `skrc` | `coulomb_rank_bounds` | [coulomb_tensor_rank_range](../../src/lib/core/coulomb_tensor_rank_range.f90) |
| `slater` | `calculate_slater_radial_integral` | [slater_radial_integral](../../src/lib/core/slater_radial_integral.f90) |
| `speak` | `buffer_mcp_coefficient` | [mcp_coefficient_output](../../src/lib/core/mcp_coefficient_output.f90) |
| `spicmv2` | `multiply_sparse_symmetric_matrix` | [sparse_symmetric_matvec](../../src/lib/core/sparse_symmetric_matvec.f90) |
| `start` | `solve_radial_origin` | [radial_origin_solution](../../src/lib/core/radial_origin_solution.f90) |
| `starttime` | `begin_run_timer` | [run_timer](../../src/lib/core/run_timer.f90) |
| `stoptime` | `finish_run_timer` | [run_timer](../../src/lib/core/run_timer.f90) |
| `yzk` | `evaluate_hartree_yz_functions` | [hartree_yz_functions](../../src/lib/core/hartree_yz_functions.f90) |

## 自旋–角向库

| 原过程名 | 当前过程名 | 源文件 |
| --- | --- | --- |
| `a1jj` | `evaluate_jj_a_few_electrons` | [jj_a_matrix_element_few_electrons](../../src/lib/spin_angular/jj_a_matrix_element_few_electrons.f90) |
| `awp1` | `evaluate_jj_aw` | [jj_aw_matrix_element](../../src/lib/spin_angular/jj_aw_matrix_element.f90) |
| `awp1jjg` | `evaluate_high_j_aw` | [jj_aw_matrix_element_high_j](../../src/lib/spin_angular/jj_aw_matrix_element_high_j.f90) |
| `c0t5s` | `evaluate_quasispin_cg_rank_half` | [quasispin_cg_rank_half](../../src/lib/spin_angular/quasispin_cg_rank_half.f90) |
| `c1e0sm` | `evaluate_quasispin_cg_rank_one_zero_projection` | [quasispin_cg_rank_one_zero_projection](../../src/lib/spin_angular/quasispin_cg_rank_one_zero_projection.f90) |
| `c1e1sm` | `evaluate_quasispin_cg_rank_one_unit_projection` | [quasispin_cg_rank_one_unit_projection](../../src/lib/spin_angular/quasispin_cg_rank_one_unit_projection.f90) |
| `cle0sm` | `evaluate_quasispin_cg_zero_projection` | [quasispin_cg_zero_projection](../../src/lib/spin_angular/quasispin_cg_zero_projection.f90) |
| `coulom` | `coulomb_angular_factor_value` | [coulomb_angular_factor](../../src/lib/spin_angular/coulomb_angular_factor.f90) |
| `diaga1` | `recouple_first_active_shell` | [first_active_shell_recoupling](../../src/lib/spin_angular/first_active_shell_recoupling.f90) |
| `diaga2` | `recouple_later_active_shell` | [later_active_shell_recoupling](../../src/lib/spin_angular/later_active_shell_recoupling.f90) |
| `diaga3` | `recouple_spectator_shell_chain` | [spectator_shell_chain_recoupling](../../src/lib/spin_angular/spectator_shell_chain_recoupling.f90) |
| `diaga4` | `recouple_pair_tensor` | [pair_tensor_recoupling](../../src/lib/spin_angular/pair_tensor_recoupling.f90) |
| `diaga5` | `recouple_final_shell` | [final_shell_recoupling](../../src/lib/spin_angular/final_shell_recoupling.f90) |
| `dracah` | `evaluate_racah_coefficient` | [racah_coefficients](../../src/lib/spin_angular/racah_coefficients.f90) |
| `eile` | `order_three_subshells` | [three_subshell_ordering](../../src/lib/spin_angular/three_subshell_ordering.f90) |
| `el1` | `accumulate_two_electron_unchanged_occupations` | [two_electron_unchanged_occupations](../../src/lib/spin_angular/two_electron_unchanged_occupations.f90) |
| `el2` | `accumulate_two_electron_pair_transfer` | [two_electron_pair_transfer](../../src/lib/spin_angular/two_electron_pair_transfer.f90) |
| `el3` | `accumulate_two_electron_single_transfer` | [two_electron_single_transfer](../../src/lib/spin_angular/two_electron_single_transfer.f90) |
| `el31` | `accumulate_single_transfer_triple_on_first_shell` | [single_transfer_triple_on_first_shell](../../src/lib/spin_angular/single_transfer_triple_on_first_shell.f90) |
| `el32` | `accumulate_single_transfer_triple_on_second_shell` | [single_transfer_triple_on_second_shell](../../src/lib/spin_angular/single_transfer_triple_on_second_shell.f90) |
| `el33` | `accumulate_single_transfer_with_spectator_shell` | [single_transfer_with_spectator_shell](../../src/lib/spin_angular/single_transfer_with_spectator_shell.f90) |
| `el4` | `accumulate_three_shell_transfers` | [two_electron_three_shell_transfer](../../src/lib/spin_angular/two_electron_three_shell_transfer.f90) |
| `el41` | `accumulate_three_shell_pair_split_or_join` | [three_shell_pair_split_or_join](../../src/lib/spin_angular/three_shell_pair_split_or_join.f90) |
| `el5` | `accumulate_four_shell_transfers` | [two_electron_four_shell_transfer](../../src/lib/spin_angular/two_electron_four_shell_transfer.f90) |
| `el51` | `accumulate_four_shell_transfer_12_34` | [four_shell_transfer_12_34](../../src/lib/spin_angular/four_shell_transfer_12_34.f90) |
| `el52` | `accumulate_four_shell_transfer_13_24` | [four_shell_transfer_13_24](../../src/lib/spin_angular/four_shell_transfer_13_24.f90) |
| `el53` | `accumulate_four_shell_transfer_14_23` | [four_shell_transfer_14_23](../../src/lib/spin_angular/four_shell_transfer_14_23.f90) |
| `gg1112` | `evaluate_subshell_product_wa_a` | [subshell_product_wa_a](../../src/lib/spin_angular/subshell_product_wa_a.f90) |
| `gg1122` | `evaluate_subshell_product_w_w` | [subshell_product_w_w](../../src/lib/spin_angular/subshell_product_w_w.f90) |
| `gg12` | `evaluate_subshell_product_a_a` | [subshell_product_a_a](../../src/lib/spin_angular/subshell_product_a_a.f90) |
| `gg1222` | `evaluate_subshell_product_a_aw` | [subshell_product_a_aw](../../src/lib/spin_angular/subshell_product_a_aw.f90) |
| `gg1233` | `evaluate_subshell_product_a_a_w` | [subshell_product_a_a_w](../../src/lib/spin_angular/subshell_product_a_a_w.f90) |
| `gg1234` | `evaluate_subshell_product_a_a_a_a` | [subshell_product_a_a_a_a](../../src/lib/spin_angular/subshell_product_a_a_a_a.f90) |
| `gracah1` | `scott_fischer_racah_coefficient` | [racah_coefficients_scott_fischer](../../src/lib/spin_angular/racah_coefficients_scott_fischer.f90) |
| `itjj` | `select_same_parity_intermediate_terms` | [jj_intermediate_terms_same_parity](../../src/lib/spin_angular/jj_intermediate_terms_same_parity.f90) |
| `itjj2` | `select_ket_parity_intermediate_terms` | [jj_intermediate_terms_ket_parity](../../src/lib/spin_angular/jj_intermediate_terms_ket_parity.f90) |
| `itjj3` | `select_bra_parity_intermediate_terms` | [jj_intermediate_terms_bra_parity](../../src/lib/spin_angular/jj_intermediate_terms_bra_parity.f90) |
| `itrexg` | `shared_tensor_rank_bounds` | [shared_tensor_rank_range](../../src/lib/spin_angular/shared_tensor_rank_range.f90) |
| `ittk` | `satisfies_angular_triangle_and_parity` | [angular_triangle_and_parity](../../src/lib/spin_angular/angular_triangle_and_parity.f90) |
| `ixjtik` | `allows_wigner_6j_triangles` | [wigner_6j_triangle_selection](../../src/lib/spin_angular/wigner_6j_triangle_selection.f90) |
| `izas1` | `allows_quasispin_projection` | [quasispin_projection_selection](../../src/lib/spin_angular/quasispin_projection_selection.f90) |
| `jfaze` | `subshell_permutation_phase` | [subshell_permutation_sign](../../src/lib/spin_angular/subshell_permutation_sign.f90) |
| `jthn` | `extract_packed_integer_digit` | [packed_integer_digit](../../src/lib/spin_angular/packed_integer_digit.f90) |
| `mes` | `report_subshell_tensor_error` | [subshell_tensor_diagnostics](../../src/lib/spin_angular/subshell_tensor_diagnostics.f90) |
| `nine` | `wigner_9j_value` | [wigner_9j](../../src/lib/spin_angular/wigner_9j.f90) |
| `nine0` | `wigner_9j_j_zero` | [wigner_9j_with_zero](../../src/lib/spin_angular/wigner_9j_with_zero.f90) |
| `nmtejj` | `lookup_jj_subshell_term` | [jj_subshell_term_lookup](../../src/lib/spin_angular/jj_subshell_term_lookup.f90) |
| `oneparticlejj` | `calculate_one_electron_tensor_coefficients` | [one_electron_tensor_coefficients](../../src/lib/spin_angular/one_electron_tensor_coefficients.f90) |
| `oneparticlejj1` | `calculate_tensor_same_subshell_coefficient` | [one_electron_tensor_same_subshell](../../src/lib/spin_angular/one_electron_tensor_same_subshell.f90) |
| `oneparticlejj2` | `calculate_tensor_transfer_coefficient` | [one_electron_tensor_between_subshells](../../src/lib/spin_angular/one_electron_tensor_between_subshells.f90) |
| `onescalar` | `calculate_one_electron_scalar_coefficients` | [one_electron_scalar_coefficients](../../src/lib/spin_angular/one_electron_scalar_coefficients.f90) |
| `onescalar1` | `calculate_scalar_same_subshell_coefficient` | [one_electron_scalar_same_subshell](../../src/lib/spin_angular/one_electron_scalar_same_subshell.f90) |
| `onescalar2` | `calculate_scalar_transfer_coefficient` | [one_electron_scalar_between_subshells](../../src/lib/spin_angular/one_electron_scalar_between_subshells.f90) |
| `perko1` | `pack_jj_subshell_state` | [jj_subshell_state_packing](../../src/lib/spin_angular/jj_subshell_state_packing.f90) |
| `perko2` | `prepare_jj_active_subshell_states` | [jj_active_subshell_states](../../src/lib/spin_angular/jj_active_subshell_states.f90) |
| `rec3` | `recouple_permuted_three_subshells` | [permuted_three_subshell_recoupling](../../src/lib/spin_angular/permuted_three_subshell_recoupling.f90) |
| `reco` | `select_scalar_interaction_recoupling` | [scalar_interaction_recoupling_selection](../../src/lib/spin_angular/scalar_interaction_recoupling_selection.f90) |
| `reco2` | `recouple_two_subshells` | [two_subshell_recoupling](../../src/lib/spin_angular/two_subshell_recoupling.f90) |
| `reco3` | `recouple_three_subshells` | [three_subshell_recoupling](../../src/lib/spin_angular/three_subshell_recoupling.f90) |
| `reco4` | `recouple_four_subshells` | [four_subshell_recoupling](../../src/lib/spin_angular/four_subshell_recoupling.f90) |
| `recoonescalar` | `select_one_electron_scalar_recoupling` | [one_electron_scalar_recoupling_selection](../../src/lib/spin_angular/one_electron_scalar_recoupling_selection.f90) |
| `recop00` | `select_one_electron_tensor_recoupling` | [one_electron_tensor_recoupling_selection](../../src/lib/spin_angular/one_electron_tensor_recoupling_selection.f90) |
| `recop1` | `recouple_one_electron_same_subshell` | [one_electron_same_subshell_recoupling](../../src/lib/spin_angular/one_electron_same_subshell_recoupling.f90) |
| `recop2` | `recouple_one_electron_two_subshells` | [one_electron_two_subshell_recoupling](../../src/lib/spin_angular/one_electron_two_subshell_recoupling.f90) |
| `rkco_gg` | `calculate_two_electron_angular_coefficients` | [two_electron_angular_coefficients](../../src/lib/spin_angular/two_electron_angular_coefficients.f90) |
| `rmeajj` | `evaluate_jj_reduced_cfp` | [jj_reduced_cfp](../../src/lib/spin_angular/jj_reduced_cfp.f90) |
| `rmeajj11` | `evaluate_two_electron_reduced_cfp` | [jj_reduced_cfp_two_electron_terms](../../src/lib/spin_angular/jj_reduced_cfp_two_electron_terms.f90) |
| `rmeajj9` | `evaluate_reduced_cfp_j_nine_halves` | [jj_reduced_cfp_j_nine_halves](../../src/lib/spin_angular/jj_reduced_cfp_j_nine_halves.f90) |
| `rmew1jj` | `evaluate_reduced_w_j_one_half` | [jj_reduced_w_j_one_half](../../src/lib/spin_angular/jj_reduced_w_j_one_half.f90) |
| `rmew3jj` | `evaluate_reduced_w_j_three_halves` | [jj_reduced_w_j_three_halves](../../src/lib/spin_angular/jj_reduced_w_j_three_halves.f90) |
| `rmew5jj` | `evaluate_reduced_w_j_five_halves` | [jj_reduced_w_j_five_halves](../../src/lib/spin_angular/jj_reduced_w_j_five_halves.f90) |
| `rmew7bjj` | `evaluate_reduced_w_j_seven_halves_high_rank` | [jj_reduced_w_j_seven_halves_higher_ranks](../../src/lib/spin_angular/jj_reduced_w_j_seven_halves_higher_ranks.f90) |
| `rmew7jj` | `evaluate_reduced_w_j_seven_halves` | [jj_reduced_w_j_seven_halves](../../src/lib/spin_angular/jj_reduced_w_j_seven_halves.f90) |
| `rumtjj` | `decode_jj_subshell_term` | [jj_subshell_term_decoding](../../src/lib/spin_angular/jj_subshell_term_decoding.f90) |
| `rwjj` | `evaluate_jj_reduced_w` | [jj_reduced_w](../../src/lib/spin_angular/jj_reduced_w.f90) |
| `sixj` | `wigner_6j_value` | [wigner_6j](../../src/lib/spin_angular/wigner_6j.f90) |
| `sixj1` | `wigner_6j_j_one` | [wigner_6j_with_j_one](../../src/lib/spin_angular/wigner_6j_with_j_one.f90) |
| `sixj2` | `wigner_6j_j_two` | [wigner_6j_with_j_two](../../src/lib/spin_angular/wigner_6j_with_j_two.f90) |
| `sixj3` | `wigner_6j_j_three` | [wigner_6j_with_j_three](../../src/lib/spin_angular/wigner_6j_with_j_three.f90) |
| `sixj35` | `wigner_6j_j_three_halves` | [wigner_6j_with_j_three_halves](../../src/lib/spin_angular/wigner_6j_with_j_three_halves.f90) |
| `sixj4` | `wigner_6j_j_four` | [wigner_6j_with_j_four](../../src/lib/spin_angular/wigner_6j_with_j_four.f90) |
| `sixj5` | `wigner_6j_j_one_half` | [wigner_6j_with_j_one_half](../../src/lib/spin_angular/wigner_6j_with_j_one_half.f90) |
| `snrc` | `breit_rank_bounds` | [breit_radial_rank_ranges](../../src/lib/spin_angular/breit_radial_rank_ranges.f90) |
| `suwjj` | `evaluate_reduced_w_from_cfp` | [jj_reduced_w_from_cfp](../../src/lib/spin_angular/jj_reduced_w_from_cfp.f90) |
| `w1jjg` | `evaluate_high_j_w` | [jj_w_matrix_element_high_j](../../src/lib/spin_angular/jj_w_matrix_element_high_j.f90) |
| `wap1` | `evaluate_jj_wa` | [jj_wa_matrix_element](../../src/lib/spin_angular/jj_wa_matrix_element.f90) |
| `wap1jjg` | `evaluate_high_j_wa` | [jj_wa_matrix_element_high_j](../../src/lib/spin_angular/jj_wa_matrix_element_high_j.f90) |
| `wj1` | `evaluate_jj_w` | [jj_w_matrix_element](../../src/lib/spin_angular/jj_w_matrix_element.f90) |
| `ww1` | `evaluate_jj_scalar_ww` | [jj_ww_scalar_matrix_element](../../src/lib/spin_angular/jj_ww_scalar_matrix_element.f90) |

## Davidson 本征求解库

| 原过程名 | 当前过程名 | 源文件 |
| --- | --- | --- |
| `adds` | `update_davidson_projected_matrix` | [davidson_projected_matrix](../../src/lib/davidson_eigensolver/davidson_projected_matrix.f90) |
| `dvdrvr` | `iterate_davidson_solver` | [davidson_iteration](../../src/lib/davidson_eigensolver/davidson_iteration.f90) |
| `dvdson` | `control_davidson_solver` | [davidson_solver_control](../../src/lib/davidson_eigensolver/davidson_solver_control.f90) |
| `gdvd` | `solve_davidson_eigenpairs` | [davidson_eigensolver](../../src/lib/davidson_eigensolver/davidson_eigensolver.f90) |
| `iniest` | `estimate_davidson_sparse_eigenpairs` | [davidson_sparse_initial_estimates](../../src/lib/davidson_eigensolver/davidson_sparse_initial_estimates.f90) |
| `initdvd` | `initialize_davidson_subspace` | [davidson_initial_subspace](../../src/lib/davidson_eigensolver/davidson_initial_subspace.f90) |
| `mgs_nrm` | `orthonormalize_davidson_vectors` | [davidson_orthonormalization](../../src/lib/davidson_eigensolver/davidson_orthonormalization.f90) |
| `multbc` | `transform_davidson_basis` | [davidson_basis_transform](../../src/lib/davidson_eigensolver/davidson_basis_transform.f90) |
| `newvec` | `expand_davidson_residuals` | [davidson_residual_expansion](../../src/lib/davidson_eigensolver/davidson_residual_expansion.f90) |
| `ovflow` | `restart_davidson_subspace` | [davidson_subspace_restart](../../src/lib/davidson_eigensolver/davidson_subspace_restart.f90) |
| `tstsel` | `check_davidson_convergence` | [davidson_convergence](../../src/lib/davidson_eigensolver/davidson_convergence.f90) |

## MPI 支持库

| 原过程名 | 当前过程名 | 源文件 |
| --- | --- | --- |
| `cpath` | `prepare_mpi_work_directories` | [mpi_work_paths](../../src/lib/mpi/mpi_work_paths.f90) |
| `gdmaxmpi_root` | `reduce_mpi_real_max_to_root` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `gdrsummpi_root` | `reduce_mpi_real_sum_to_root` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `gdsummpi` | `allreduce_mpi_real_sum` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `gisummpi` | `allreduce_mpi_integer_sum` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `iniestmpi` | `estimate_mpi_initial_eigenvectors` | [mpi_initial_eigenvectors](../../src/lib/mpi/mpi_initial_eigenvectors.f90) |
| `lodrwfmpi` | `read_mpi_radial_wavefunction_records` | [mpi_radial_wavefunction_reader](../../src/lib/mpi/mpi_radial_wavefunction_reader.f90) |
| `mpix_chkpt` | `check_mpi_checkpoint` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `mpix_printmsg` | `print_mpi_message` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `setisompi` | `read_mpi_isotope_data` | [mpi_isotope_input](../../src/lib/mpi/mpi_isotope_input.f90) |
| `setrwfmpi` | `read_mpi_radial_wavefunction_file` | [mpi_radial_wavefunction_file](../../src/lib/mpi/mpi_radial_wavefunction_file.f90) |
| `spicmvmpi` | `multiply_mpi_sparse_matrix` | [mpi_sparse_matrix_product](../../src/lib/mpi/mpi_sparse_matrix_product.f90) |
| `startmpi` | `initialize_mpi_runtime` | [mpi_runtime](../../src/lib/mpi/mpi_runtime.f90) |
| `startmpi2` | `initialize_mpi_runtime_with_status` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `stopmpi` | `finalize_mpi_runtime` | [mpi_runtime](../../src/lib/mpi/mpi_runtime.f90) |
| `stopmpi2` | `finish_mpi_application` | [mpi_utilities](../../src/lib/mpi/mpi_utilities.f90) |
| `sys_chdir` | `change_mpi_directory` | [mpi_change_directory](../../src/lib/mpi/mpi_change_directory.f90) |
| `sys_getwd` | `get_mpi_working_directory` | [mpi_current_directory](../../src/lib/mpi/mpi_current_directory.f90) |
| `sys_mkdir` | `create_mpi_directory` | [mpi_make_directory](../../src/lib/mpi/mpi_make_directory.f90) |

## Breit 角向支持库

`breit_angular_coefficient_buffer` 从 `store_breit_angular_coefficient` 派生外部
`collect_breit_angular_coefficient` 的显式接口。库默认实现和 RCI 替换实现都使用新链接符号。
原先调用或定义旧入口的 Fortran 客户端需更新过程名，重新编译并链接。

| 原过程名 | 当前过程名 | 源文件 |
| --- | --- | --- |
| `talk`、`buffer_breit_angular_coefficient` | `collect_breit_angular_coefficient` | [breit_angular_coefficient_buffer](../../src/lib/breit_angular_support/breit_angular_coefficient_buffer.f90) |
| `talk_fallback` | `store_breit_angular_coefficient` | [breit_angular_coefficient_buffer](../../src/lib/breit_angular_support/breit_angular_coefficient_buffer.f90) |
