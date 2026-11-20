# GRASP 公共基础库命名迁移

日期 2026-11-20。源码基准为 `2d4c721` 加同一工作树中的变量和接口整理。

| 用途 | 原名 | 当前名称 |
| --- | --- | --- |
| 源码目录 | `src/lib/lib9290` | [src/lib/core](../../src/lib/core) |
| CMake 目标 | `9290` | `core` |
| 静态库 | `lib9290.a` | `libcore.a` |
| Fortran 模块目录（构建树） | `modules/9290` | `modules/core` |
| 径向积分回归 | `lib9290_quad` | `core_quadrature` |

源码文件名与主模块名一致，使用带下划线的小写职责名称，例如 `radial_grid_setup` 和 `csf_quantum_queries`。在套件外使用这些模块的 Fortran 代码，需更新 `USE` 导入，然后重新编译并链接 `libcore.a`。CMake 调用方可使用 `core` 目标，继续由它传递 Rust CSF 描述符库依赖。

GNU Fortran 15 的干净 Release 构建启用 MPI 和 OpenMP，110 项 CTest 全部通过。关闭 MPI 和 OpenMP 的 Check 构建通过 70 项测试，Rust 读库测试在 Release 中执行。安装验证使用 `libcore.a` 和安装目录中的 `csf_quantum_queries.mod` 编译、链接并运行独立调用程序。

过程名见[共享库过程对照](shared_library_procedure_names.md)。

## 文件与模块对照

| 原文件 | 原模块 | 现模块与源码 |
| --- | --- | --- |
| `alcbuf.f90` | `alcbuf_i` | `angular_coefficient_buffer` · [angular_coefficient_buffer.f90](../../src/lib/core/angular_coefficient_buffer.f90) |
| `arctan.f90` | `arctan_i` | `nonnegative_atan2` · [nonnegative_atan2.f90](../../src/lib/core/nonnegative_atan2.f90) |
| `calen.f90` | `calen_i` | `calendar_clock` · [calendar_clock.f90](../../src/lib/core/calendar_clock.f90) |
| `cgamma.f90` | `cgamma_i` | `complex_gamma` · [complex_gamma.f90](../../src/lib/core/complex_gamma.f90) |
| `clrx.f90` | `clrx_i` | `relativistic_three_j` · [relativistic_three_j.f90](../../src/lib/core/relativistic_three_j.f90) |
| `convrt.f90` | `convrt_i` | `integer_format` · [integer_format.f90](../../src/lib/core/integer_format.f90) |
| `convrt2.f90` | `convrt2_i` | `integer_format_with_context` · [integer_format_with_context.f90](../../src/lib/core/integer_format_with_context.f90) |
| `convrt_double.f90` | `convrt_double_i` | `half_integer_format` · [half_integer_format.f90](../../src/lib/core/half_integer_format.f90) |
| `cord.f90` | `cord_i` | `closed_shell_mcp` · [closed_shell_mcp.f90](../../src/lib/core/closed_shell_mcp.f90) |
| `count.f90` | `count_implementation`、`count_i` | `radial_node_count` · [radial_node_count.f90](../../src/lib/core/radial_node_count.f90) |
| `cre.f90` | `cre_i` | `reduced_angular_matrix_element` · [reduced_angular_matrix_element.f90](../../src/lib/core/reduced_angular_matrix_element.f90) |
| `csf_basis_input.f90` | `csf_basis_input` | `csf_full_basis_reader` · [csf_full_basis_reader.f90](../../src/lib/core/csf_full_basis_reader.f90) |
| `csf_basis_merge.f90` | `csf_basis_merge` | `csf_basis_merging` · [csf_basis_merging.f90](../../src/lib/core/csf_basis_merging.f90) |
| `csf_block_input.f90` | `csf_block_input` | `rci_csf_block_reader` · [rci_csf_block_reader.f90](../../src/lib/core/rci_csf_block_reader.f90) |
| `csf_descriptor_input.f90` | `csf_descriptor_input` | `csf_descriptor_reader` · [csf_descriptor_reader.f90](../../src/lib/core/csf_descriptor_reader.f90) |
| `csf_queries.f90` | `csf_queries` | `csf_quantum_queries` · [csf_quantum_queries.f90](../../src/lib/core/csf_quantum_queries.f90) |
| `csf_rust_ffi.f90` | `csf_rust_ffi` | `csf_descriptor_c_bindings` · [csf_descriptor_c_bindings.f90](../../src/lib/core/csf_descriptor_c_bindings.f90) |
| `dcbsrw.f90` | `dcbsrw_i` | `dirac_coulomb_radial_orbital` · [dirac_coulomb_radial_orbital.f90](../../src/lib/core/dirac_coulomb_radial_orbital.f90) |
| `dinit.f90` | `dinit_i` | `real_vector_initialization` · [real_vector_initialization.f90](../../src/lib/core/real_vector_initialization.f90) |
| `dmerge.f90` | `dmerge_i` | `sparse_dot_and_scatter` · [sparse_dot_and_scatter.f90](../../src/lib/core/sparse_dot_and_scatter.f90) |
| `dpbdt.f90` | `dpbdt_i` | `radial_wavefunction_derivatives` · [radial_wavefunction_derivatives.f90](../../src/lib/core/radial_wavefunction_derivatives.f90) |
| `draw.f90` | `draw_i` | `radial_ascii_plot` · [radial_ascii_plot.f90](../../src/lib/core/radial_ascii_plot.f90) |
| `es.f90` | `es_i` | `alternating_exponential_series` · [alternating_exponential_series.f90](../../src/lib/core/alternating_exponential_series.f90) |
| `factt.f90` | `factt_i` | `log_factorial_table` · [log_factorial_table.f90](../../src/lib/core/log_factorial_table.f90) |
| `getrsl.f90` | `getrsl_i` | `subshell_label_selection` · [subshell_label_selection.f90](../../src/lib/core/subshell_label_selection.f90) |
| `getyn.f90` | `getyn_i` | `yes_no_input` · [yes_no_input.f90](../../src/lib/core/yes_no_input.f90) |
| `ichkq1.f90` | `ichkq1_i` | `one_electron_csf_filter` · [one_electron_csf_filter.f90](../../src/lib/core/one_electron_csf_filter.f90) |
| `ichkq2.f90` | `ichkq2_i` | `two_electron_csf_filter` · [two_electron_csf_filter.f90](../../src/lib/core/two_electron_csf_filter.f90) |
| `ichop.f90` | `ichop_i` | `subshell_occupation_class` · [subshell_occupation_class.f90](../../src/lib/core/subshell_occupation_class.f90) |
| `icopy.f90` | `icopy_i` | `integer_vector_copy` · [integer_vector_copy.f90](../../src/lib/core/integer_vector_copy.f90) |
| `iniest2.f90` | `iniest2_i` | `sparse_eigensolver_initial_basis` · [sparse_eigensolver_initial_basis.f90](../../src/lib/core/sparse_eigensolver_initial_basis.f90) |
| `interp.f90` | `interp_i` | `aitken_interpolation` · [aitken_interpolation.f90](../../src/lib/core/aitken_interpolation.f90) |
| `intrpq.f90` | `intrpq_i` | `radial_orbital_interpolation` · [radial_orbital_interpolation.f90](../../src/lib/core/radial_orbital_interpolation.f90) |
| `items.f90` | `items_i` | `level_selection_parser` · [level_selection_parser.f90](../../src/lib/core/level_selection_parser.f90) |
| `itrig.f90` | `itrig_i` | `angular_momentum_triangle` · [angular_momentum_triangle.f90](../../src/lib/core/angular_momentum_triangle.f90) |
| `ldigit.f90` | `ldigit_i` | `decimal_digit_check` · [decimal_digit_check.f90](../../src/lib/core/decimal_digit_check.f90) |
| `lodiso.f90` | `lodiso_i` | `isotope_data_records` · [isotope_data_records.f90](../../src/lib/core/isotope_data_records.f90) |
| `lodrwf.f90` | `lodrwf_i` | `radial_wavefunction_records` · [radial_wavefunction_records.f90](../../src/lib/core/radial_wavefunction_records.f90) |
| `lodstate.f90` | `lodstate_i` | `block_level_selection` · [block_level_selection.f90](../../src/lib/core/block_level_selection.f90) |
| `ltab.f90` | `ltab_i` | `subshell_parent_table` · [subshell_parent_table.f90](../../src/lib/core/subshell_parent_table.f90) |
| `npol_init.f90` | `npol_initialization` | `nuclear_polarization_initialization` · [nuclear_polarization_initialization.f90](../../src/lib/core/nuclear_polarization_initialization.f90) |
| `npol_input.f90` | `npol_input` | `nuclear_polarization_input` · [nuclear_polarization_input.f90](../../src/lib/core/nuclear_polarization_input.f90) |
| `npol_kernel.f90` | `npol_kernel` | `nuclear_polarization_potential` · [nuclear_polarization_potential.f90](../../src/lib/core/nuclear_polarization_potential.f90) |
| `nucpot.f90` | `nucpot_i` | `nuclear_coulomb_potential` · [nuclear_coulomb_potential.f90](../../src/lib/core/nuclear_coulomb_potential.f90) |
| `openfl.f90` | `openfl_i` | `file_open` · [file_open.f90](../../src/lib/core/file_open.f90) |
| `orthsc.f90` | `orthsc_i` | `radial_orbital_orthogonalization` · [radial_orbital_orthogonalization.f90](../../src/lib/core/radial_orbital_orthogonalization.f90) |
| `pack.f90` | `pack_implementation`、`pack_i` | `csf_byte_packing` · [csf_byte_packing.f90](../../src/lib/core/csf_byte_packing.f90) |
| `posfile.f90` | `posfile_i` | `file_positioning` · [file_positioning.f90](../../src/lib/core/file_positioning.f90) |
| `quad.f90` | `quad_i` | `radial_quadrature` · [radial_quadrature.f90](../../src/lib/core/radial_quadrature.f90) |
| `radgrd.f90` | `radgrd_i` | `radial_grid_setup` · [radial_grid_setup.f90](../../src/lib/core/radial_grid_setup.f90) |
| `rint.f90` | `rint_i` | `weighted_radial_overlap` · [weighted_radial_overlap.f90](../../src/lib/core/weighted_radial_overlap.f90) |
| `rinti.f90` | `rinti_i` | `one_electron_radial_integral` · [one_electron_radial_integral.f90](../../src/lib/core/one_electron_radial_integral.f90) |
| `rmcdhf_descriptor_input.f90` | `rmcdhf_descriptor_input` | `rmcdhf_csf_input` · [rmcdhf_csf_input.f90](../../src/lib/core/rmcdhf_csf_input.f90) |
| `setcon.f90` | `setcon_i` | `physical_constants` · [physical_constants.f90](../../src/lib/core/physical_constants.f90) |
| `setcsla.f90` | `setcsla_i` | `csf_basis_loader` · [csf_basis_loader.f90](../../src/lib/core/csf_basis_loader.f90) |
| `setiso.f90` | `setiso_i` | `isotope_data_input` · [isotope_data_input.f90](../../src/lib/core/isotope_data_input.f90) |
| `setj.f90` | `setj_i` | `angular_recoupling_tables` · [angular_recoupling_tables.f90](../../src/lib/core/angular_recoupling_tables.f90) |
| `setmc.f90` | `setmc_i` | `machine_constants` · [machine_constants.f90](../../src/lib/core/machine_constants.f90) |
| `setpot.f90` | `setpot_i` | `radial_dirac_coefficients` · [radial_dirac_coefficients.f90](../../src/lib/core/radial_dirac_coefficients.f90) |
| `setqic.f90` | `setqic_i` | `radial_stencil_weights` · [radial_stencil_weights.f90](../../src/lib/core/radial_stencil_weights.f90) |
| `setqna.f90` | `setqna_i` | `csf_pair_quantum_numbers` · [csf_pair_quantum_numbers.f90](../../src/lib/core/csf_pair_quantum_numbers.f90) |
| `setrwfa.f90` | `setrwfa_i` | `radial_wavefunction_input` · [radial_wavefunction_input.f90](../../src/lib/core/radial_wavefunction_input.f90) |
| `skrc.f90` | `skrc_i` | `coulomb_tensor_rank_range` · [coulomb_tensor_rank_range.f90](../../src/lib/core/coulomb_tensor_rank_range.f90) |
| `slater.f90` | `slater_i` | `slater_radial_integral` · [slater_radial_integral.f90](../../src/lib/core/slater_radial_integral.f90) |
| `speak.f90` | `speak_i` | `mcp_coefficient_output` · [mcp_coefficient_output.f90](../../src/lib/core/mcp_coefficient_output.f90) |
| `spicmv2.f90` | `spicmv2_i` | `sparse_symmetric_matvec` · [sparse_symmetric_matvec.f90](../../src/lib/core/sparse_symmetric_matvec.f90) |
| `start.f90` | `start_i` | `radial_origin_solution` · [radial_origin_solution.f90](../../src/lib/core/radial_origin_solution.f90) |
| `starttime.f90` | `starttime_i` | `run_timer_start` · [run_timer_start.f90](../../src/lib/core/run_timer_start.f90) |
| `stoptime.f90` | `stoptime_i` | `run_timer_finish` · [run_timer_finish.f90](../../src/lib/core/run_timer_finish.f90) |
| `suite_system.f90` | `suite_system` | `posix_services` · [posix_services.f90](../../src/lib/core/posix_services.f90) |
| `yzk.f90` | `yzk_i` | `hartree_yz_functions` · [hartree_yz_functions.f90](../../src/lib/core/hartree_yz_functions.f90) |
