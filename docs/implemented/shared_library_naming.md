# GRASP 共享库命名

命名迁移基于 2026-11-20 的 suite 工作树，Git 基点为 `2d4c721`。

| 项目 | 原名 | 现名 |
| --- | --- | --- |
| 源码目录 | `src/lib/libmod/` | [`src/lib/grasp_shared/`](../../src/lib/grasp_shared/) |
| CMake 目标 | `mod` | `grasp_shared` |
| 静态库 | `libmod.a` | `libgrasp_shared.a` |
| 链接参数 | `-lmod` | `-lgrasp_shared` |
| 构建模块目录 | `<build>/modules/mod/` | `<build>/modules/grasp_shared/` |
| 安装模块目录 | `<prefix>/lib/mod/` | `<prefix>/lib/grasp_shared/` |

每个 `.f90` 文件定义一个同名模块。文件名和模块名使用小写 snake_case，按数据或功能命名。
`*_state` 表示计算过程维护的状态，`*_workspace` 表示可反复使用的工作区，`*_table` 表示按索引访问的数据。
保留 CSF、SCF、MCP、RCI、QED 和 jj/LSJ 等领域缩写。

安装路径中的 `<prefix>` 指 CMake 配置时设置的 `CMAKE_INSTALL_PREFIX`。

库名和模块名是链接接口的一部分。外部 Fortran 调用方需要按下表更新 `USE`，
使用新的 `.mod` 文件重新编译，并链接 `libgrasp_shared.a`。原模块中的数据和过程接口继续使用原有名称。
`grasp_parameters` 转导 `cal_constants` 中的六个语义容量常量。
`npar_C.f90` 中的两个独立模块各自迁移到同名文件。
备用的内存分析源码沿用新模块名，文件为 `array_memory_management.f90_PROFILING`。

GNU Fortran 15.3 构建通过。Release 启用串行、MPI 和 jj2lsj OpenMP，110 项 CTest 全部通过。
串行 Check 的 70 项测试全部通过，Rust 描述符读库测试在 Release 中执行。
安装验证核对了 109 份 `.mod`，并使用安装产物编译、链接和运行了一个独立的 Fortran 调用程序。

## 文件与模块对照

| 原文件 | 原模块 | 新模块和源文件 | 职责 |
| --- | --- | --- | --- |
| `AME_C.f90` | `ame_c` | [`angular_matrix_element_workspace`](../../src/lib/grasp_shared/angular_matrix_element_workspace.f90) | Angular matrix element indices and work arrays |
| `bcore_C.f90` | `bcore_c` | [`breit_core_occupations`](../../src/lib/grasp_shared/breit_core_occupations.f90) | Core occupations used for Breit contributions |
| `bess_C.f90` | `bess_c` | [`bessel_function_workspace`](../../src/lib/grasp_shared/bessel_function_workspace.f90) | Radial Bessel and Neumann function arrays |
| `bilst_C.f90` | `bilst_c` | [`breit_integral_tables`](../../src/lib/grasp_shared/breit_integral_tables.f90) | Indices and values for the six Breit integral types |
| `biorb_C.f90` | `biorb_c` | [`biorthogonal_orbital_sets`](../../src/lib/grasp_shared/biorthogonal_orbital_sets.f90) | Initial and final orbital and basis metadata |
| `blim_C.f90` | `blim_c` | [`rci_prerun_settings`](../../src/lib/grasp_shared/rci_prerun_settings.f90) | CSF limits and coefficient thresholds for the RCI prerun |
| `blk_C.f90` | `blk_c` | [`transition_basis_blocks`](../../src/lib/grasp_shared/transition_basis_blocks.f90) | Initial and final basis block sizes and indices for transitions |
| `blkidx_C.f90` | `blkidx_c` | [`transition_block_indices`](../../src/lib/grasp_shared/transition_block_indices.f90) | Block index array for transitions |
| `buffer_C.f90` | `buffer_c` | [`angular_coefficient_buffer_state`](../../src/lib/grasp_shared/angular_coefficient_buffer_state.f90) | Angular coefficient buffer capacity, labels, and values |
| `cffmat_C.f90` | `cffmat_c` | [`final_orbital_rotation_matrix`](../../src/lib/grasp_shared/final_orbital_rotation_matrix.f90) | Final-state orbital rotation coefficients |
| `ciimat_C.f90` | `ciimat_c` | [`initial_orbital_rotation_matrix`](../../src/lib/grasp_shared/initial_orbital_rotation_matrix.f90) | Initial-state orbital rotation coefficients |
| `cimat_C.f90` | `cimat_c` | [`biorthogonal_orbital_rotations`](../../src/lib/grasp_shared/biorthogonal_orbital_rotations.f90) | Initial and final orbital rotation coefficients |
| `cnc_C.f90` | `cnc_c` | [`radial_difference_coefficients`](../../src/lib/grasp_shared/radial_difference_coefficients.f90) | Radial finite-difference coefficients |
| `coeils_C.f90` | `coeils_c` | [`coulomb_integral_table`](../../src/lib/grasp_shared/coulomb_integral_table.f90) | Indices, values, and offsets for Coulomb integrals |
| `cons_C.f90` | `cons_c` | [`numeric_constants`](../../src/lib/grasp_shared/numeric_constants.f90) | Shared numerical constants and tolerances |
| `core_C.f90` | `core_c` | [`closed_core_subshells`](../../src/lib/grasp_shared/closed_core_subshells.f90) | Count of closed-core subshells |
| `corre_C.f90` | `corre_c` | [`correlation_orbital_flags`](../../src/lib/grasp_shared/correlation_orbital_flags.f90) | Flags identifying correlation orbitals |
| `coun_C.f90` | `coun_c` | [`radial_node_count_settings`](../../src/lib/grasp_shared/radial_node_count_settings.f90) | Threshold and diagnostic context for radial node counting |
| `couple_C.f90` | `couple_c` | [`angular_coupling_workspace`](../../src/lib/grasp_shared/angular_coupling_workspace.f90) | Triads and angular momenta for recoupling |
| `cteilsrk_C.f90` | `cteilsrk_c` | [`distributed_slater_integrals`](../../src/lib/grasp_shared/distributed_slater_integrals.f90) | Distributed Slater integral cache and local entry count |
| `cuto_C.f90` | `cuto_c` | [`transition_output_cutoff`](../../src/lib/grasp_shared/transition_output_cutoff.f90) | Cutoff for printing transition data |
| `damp_C.f90` | `damp_c` | [`scf_damping_factors`](../../src/lib/grasp_shared/scf_damping_factors.f90) | Orbital and mixing coefficient damping factors |
| `debug_C.f90` | `debug_c` | [`debug_output_settings`](../../src/lib/grasp_shared/debug_output_settings.f90) | Debug output switches and cutoff |
| `decide_C.f90` | `decide_c` | [`hamiltonian_correction_flags`](../../src/lib/grasp_shared/hamiltonian_correction_flags.f90) | Switches for Breit, QED, and mass-shift terms |
| `def_C.f90` | `def_c` | [`atomic_calculation_state`](../../src/lib/grasp_shared/atomic_calculation_state.f90) | Atomic system, unit conversions, and calculation parameters |
| `default_C.f90` | `default_c` | [`calculation_input_defaults`](../../src/lib/grasp_shared/calculation_input_defaults.f90) | Calculation defaults and dump/input flags |
| `dumx_C.f90` | `dumx_c` | [`angular_recoupling_backup`](../../src/lib/grasp_shared/angular_recoupling_backup.f90) | Saved peel-orbital indices and bra/ket coupling weights |
| `eigv_C.f90` | `eigv_c` | [`eigenstate_expansions`](../../src/lib/grasp_shared/eigenstate_expansions.f90) | Energy shifts, eigenvalues, and CSF mixing coefficients |
| `eigvec1_C.f90` | `eigvec1_c` | [`scf_previous_eigenvectors`](../../src/lib/grasp_shared/scf_previous_eigenvectors.f90) | Previous eigenvectors used by SCF |
| `facts_C.f90` | `facts_c` | [`log_factorial_cache`](../../src/lib/grasp_shared/log_factorial_cache.f90) | Shared log-factorial table |
| `fixd_C.f90` | `fixd_c` | [`fixed_orbital_flags`](../../src/lib/grasp_shared/fixed_orbital_flags.f90) | Count and flags of fixed orbitals |
| `foparm_C.f90` | `foparm_c` | [`first_order_csf_cutoffs`](../../src/lib/grasp_shared/first_order_csf_cutoffs.f90) | CSF cutoffs for first-order restrictions |
| `fposition_C.f90` | `fposition_c` | [`rci_matrix_file_position`](../../src/lib/grasp_shared/rci_matrix_file_position.f90) | Position in the RCI matrix file |
| `grid_C.f90` | `grid_c` | [`radial_grid_state`](../../src/lib/grasp_shared/radial_grid_state.f90) | Radial coordinates, derivatives, scale, steps, and point count |
| `hblock_C.f90` | `hblock_c` | [`hamiltonian_symmetry_blocks`](../../src/lib/grasp_shared/hamiltonian_symmetry_blocks.f90) | CSF and eigenvector counts by symmetry block |
| `hmat_C.f90` | `hmat_c` | [`sparse_hamiltonian_storage`](../../src/lib/grasp_shared/sparse_hamiltonian_storage.f90) | Sparse Hamiltonian entries, indices, and matrix counters |
| `horb_C.f90` | `horb_c` | [`hydrogenic_radial_components`](../../src/lib/grasp_shared/hydrogenic_radial_components.f90) | Hydrogenic large and small radial components for QED |
| `hydpar_C.f90` | `hydpar_c` | [`hydrogenic_screening`](../../src/lib/grasp_shared/hydrogenic_screening.f90) | Orbital screening parameters for hydrogenic estimates |
| `iccu_C.f90` | `iccu_c` | [`csf_block_cutoffs`](../../src/lib/grasp_shared/csf_block_cutoffs.f90) | CSF cutoff arrays for symmetry blocks |
| `input_io.f90` | `input_io` | [`input_error_handling`](../../src/lib/grasp_shared/input_error_handling.f90) | Input error reporting and termination callback |
| `int_C.f90` | `int_c` | [`radial_dirac_workspace`](../../src/lib/grasp_shared/radial_dirac_workspace.f90) | Working radial orbital, boundary data, and solver arrays |
| `invt_C.f90` | `invt_c` | [`orbital_inversion_flags`](../../src/lib/grasp_shared/orbital_inversion_flags.f90) | Per-orbital flags for inversion of the radial solution |
| `iounit_C.f90` | `iounit_c` | [`standard_io_units`](../../src/lib/grasp_shared/standard_io_units.f90) | Standard input, output, and error unit numbers |
| `jj2lsj_C.f90` | `jj2lsj_c` | [`ls_jj_subshell_terms`](../../src/lib/grasp_shared/ls_jj_subshell_terms.f90) | LS subshell terms and LS–jj transformation data types |
| `jj2lsjbio_C.f90` | `jj2lsjbio_c` | [`transition_lsj_labels`](../../src/lib/grasp_shared/transition_lsj_labels.f90) | LSJ level labels, energies, and indices for transitions |
| `jlabl_C.f90` | `jlabl_c` | [`angular_symmetry_labels`](../../src/lib/grasp_shared/angular_symmetry_labels.f90) | Text labels for angular momentum and parity |
| `jqjc_C.f90` | `jqjc_c` | [`biorthogonal_transition_selection`](../../src/lib/grasp_shared/biorthogonal_transition_selection.f90) | Initial and final angular symmetry selections |
| `keilst_C.f90` | `keilst_c` | [`kinetic_energy_integral_table`](../../src/lib/grasp_shared/kinetic_energy_integral_table.f90) | Indices, values, and offsets for kinetic-energy integrals |
| `kkstart_C.f90` | `kkstart_c` | [`coulomb_rank_offsets`](../../src/lib/grasp_shared/coulomb_rank_offsets.f90) | Start offsets by Coulomb multipole rank |
| `kkstartbreit_C.f90` | `kkstartbreit_c` | [`breit_rank_offsets`](../../src/lib/grasp_shared/breit_rank_offsets.f90) | Start offsets by rank for the two Breit groups |
| `l1_C.f90` | `l1_c` | [`csf_pair_subshell_couplings`](../../src/lib/grasp_shared/csf_pair_subshell_couplings.f90) | Bra/ket subshell and total coupling quantum numbers |
| `l2_C.f90` | `l2_c` | [`angular_coupling_triads`](../../src/lib/grasp_shared/angular_coupling_triads.f90) | Stored recoupling triads |
| `lagr_C.f90` | `lagr_c` | [`orbital_lagrange_multipliers`](../../src/lib/grasp_shared/orbital_lagrange_multipliers.f90) | Lagrange multipliers and their orbital indices |
| `left_C.f90` | `left_c` | [`initial_orbital_selection`](../../src/lib/grasp_shared/initial_orbital_selection.f90) | Selection flags for constructing initial orbitals |
| `lib92p_C.f90` | `lib92p_c` | [`legacy_plant_settings`](../../src/lib/grasp_shared/legacy_plant_settings.f90) | Legacy NPLANT and LPLANT settings |
| `lic13_C.f90` | `lic13_c` | [`radial_derivative_stencil`](../../src/lib/grasp_shared/radial_derivative_stencil.f90) | Matrix of coefficients for the radial derivative stencil |
| `linear_algebra_interfaces.f90` | `suite_linear_algebra` | [`blas_lapack_interfaces`](../../src/lib/grasp_shared/blas_lapack_interfaces.f90) | Explicit interfaces for external BLAS and LAPACK routines |
| `m_C.f90` | `m_c` | [`csf_pair_angular_state`](../../src/lib/grasp_shared/csf_pair_angular_state.f90) | Bra/ket occupations and couplings plus peel/core orbital lists |
| `mcp_C.f90` | `mcp_c` | [`mcp_generation_options`](../../src/lib/grasp_shared/mcp_generation_options.f90) | Maximum Coulomb multipole rank and MCP generation options |
| `mcpa_C.f90` | `mcpa_c` | [`mcp_file_metadata`](../../src/lib/grasp_shared/mcp_file_metadata.f90) | Maximum multipole rank in the MCP files |
| `mcpb_C.f90` | `mcpb_c` | [`rmcdhf_mcp_options`](../../src/lib/grasp_shared/rmcdhf_mcp_options.f90) | MCP diagonal and first-order settings for RMCDHF |
| `mcpdata_C.f90` | `mcpdata_c` | [`mcp_coefficient_storage`](../../src/lib/grasp_shared/mcp_coefficient_storage.f90) | MCP coefficients, CSF pairs, and integral indices |
| `memory_man.f90` | `memory_man` | [`array_memory_management`](../../src/lib/grasp_shared/array_memory_management.f90) | Allocation, resizing, and deallocation wrappers |
| `mtjj2_C.f90` | `mtjj2_c` | [`high_j_subshell_term_data`](../../src/lib/grasp_shared/high_j_subshell_term_data.f90) | Subshell term data for twice-j values 9 and 11 |
| `mtjj_C.f90` | `mtjj_c` | [`subshell_term_multiplicities`](../../src/lib/grasp_shared/subshell_term_multiplicities.f90) | Multiplicity data for jj-coupled subshell terms |
| `ncc_C.f90` | `ncc_c` | [`radial_quadrature_weights`](../../src/lib/grasp_shared/radial_quadrature_weights.f90) | Weights used by radial quadrature and related stencils |
| `ncdist_C.f90` | `ncdist_c` | [`nuclear_charge_distribution`](../../src/lib/grasp_shared/nuclear_charge_distribution.f90) | Radial nuclear charge distribution array |
| `node_C.f90` | `node_c` | [`expected_orbital_nodes`](../../src/lib/grasp_shared/expected_orbital_nodes.f90) | Expected radial node count per orbital |
| `npar_C.f90` | `npar_c` | [`nuclear_model_parameters_state`](../../src/lib/grasp_shared/nuclear_model_parameters_state.f90) | Parameters of the nuclear charge model |
| `npar_C.f90` | `nsmdat_c` | [`isotope_properties`](../../src/lib/grasp_shared/isotope_properties.f90) | Nuclear spin, moments, and isotope property factors |
| `npol_C.f90` | `npol_c` | [`nuclear_polarization_state`](../../src/lib/grasp_shared/nuclear_polarization_state.f90) | Nuclear polarization configuration and published radial potentials |
| `npot_C.f90` | `npot_c` | [`nuclear_potential_state`](../../src/lib/grasp_shared/nuclear_potential_state.f90) | Nuclear Coulomb potential array and Fermi radius index |
| `offd_C.f90` | `offd_c` | [`transition_off_diagonal_flags`](../../src/lib/grasp_shared/transition_off_diagonal_flags.f90) | Off-diagonal transition calculation switches |
| `orb_C.f90` | `orb_c` | [`orbital_csf_basis`](../../src/lib/grasp_shared/orbital_csf_basis.f90) | Orbital quantum numbers, energies, and CSF occupations |
| `orba_C.f90` | `orba_c` | [`scf_orbital_order`](../../src/lib/grasp_shared/scf_orbital_order.f90) | Orbital processing order during SCF |
| `orbord_C.f90` | `orbord_c` | [`biorthogonal_orbital_order`](../../src/lib/grasp_shared/biorthogonal_orbital_order.f90) | Initial and final orbital ordering flags |
| `orthct_C.f90` | `orthct_c` | [`scf_orthogonality_control`](../../src/lib/grasp_shared/scf_orthogonality_control.f90) | SCF orthogonality enforcement flag |
| `osc_C.f90` | `osc_c` | [`transition_matrix_workspace`](../../src/lib/grasp_shared/transition_matrix_workspace.f90) | Angular, radial, and matrix arrays for transition calculations |
| `ovl_C.f90` | `ovl_c` | [`orbital_overlap_control`](../../src/lib/grasp_shared/orbital_overlap_control.f90) | Orbital overlap control flag |
| `parameter_def_M.f90` | `parameter_def` | [`grasp_parameters`](../../src/lib/grasp_shared/grasp_parameters.f90) | Public re-export of the six centrally defined capacity parameters |
| `peav_C.f90` | `peav_c` | [`block_mean_energies`](../../src/lib/grasp_shared/block_mean_energies.f90) | Mean energy array by symmetry block |
| `pos_C.f90` | `pos_c` | [`symmetry_block_offsets`](../../src/lib/grasp_shared/symmetry_block_offsets.f90) | Offsets into CSF, eigenvalue, and mixing coefficient arrays |
| `pote_C.f90` | `pote_c` | [`scf_radial_potentials`](../../src/lib/grasp_shared/scf_radial_potentials.f90) | Direct and exchange radial potential arrays |
| `prnt_C.f90` | `prnt_c` | [`selected_level_indices`](../../src/lib/grasp_shared/selected_level_indices.f90) | Selected eigenvectors and indices for the active and paired bases |
| `qedcut_C.f90` | `qedcut_c` | [`qed_orbital_cutoffs`](../../src/lib/grasp_shared/qed_orbital_cutoffs.f90) | Orbital cutoffs used for QED corrections |
| `rang_Int_C.f90` | `rang_int_c` | [`csf_match_workspace`](../../src/lib/grasp_shared/csf_match_workspace.f90) | CSF match results and shell and coupling strings |
| `ribojj11_C.f90` | `ribojj11_c` | [`jj_reduced_coefficients_11`](../../src/lib/grasp_shared/jj_reduced_coefficients_11.f90) | Reduced-coefficient tables for twice-j equal to 11 |
| `ribojj9_C.f90` | `ribojj9_c` | [`jj_reduced_coefficients_9`](../../src/lib/grasp_shared/jj_reduced_coefficients_9.f90) | Reduced-coefficient tables for twice-j equal to 9 |
| `ribojj_C.f90` | `ribojj_c` | [`jj_reduced_coefficients`](../../src/lib/grasp_shared/jj_reduced_coefficients.f90) | Reduced-coefficient tables for jj-coupled subshells |
| `sacoef_C.f90` | `sacoef_c` | [`angular_coefficient_channels`](../../src/lib/grasp_shared/angular_coefficient_channels.f90) | Angular coefficients and CSF and orbital labels by channel |
| `sbc_C.f90` | `sbc_c` | [`radial_stepping_coefficients`](../../src/lib/grasp_shared/radial_stepping_coefficients.f90) | Coefficients for radial integration steps |
| `sbdat1_C.f90` | `sbdat1_c` | [`biorthogonal_subshell_mapping`](../../src/lib/grasp_shared/biorthogonal_subshell_mapping.f90) | Subshell lookup arrays for the biorthogonal transformation |
| `sbdat_C.f90` | `sbdat_c` | [`biorthogonal_subshell_layout`](../../src/lib/grasp_shared/biorthogonal_subshell_layout.f90) | Initial and final subshell indices grouped by symmetry |
| `scf_C.f90` | `scf_c` | [`scf_iteration_state`](../../src/lib/grasp_shared/scf_iteration_state.f90) | SCF solver controls and coefficient work arrays |
| `stat_C.f90` | `stat_c` | [`csf_coupling_data`](../../src/lib/grasp_shared/csf_coupling_data.f90) | Stored subshell and total CSF coupling codes |
| `stor_C.f90` | `stor_c` | [`breit_bessel_cache_keys`](../../src/lib/grasp_shared/breit_bessel_cache_keys.f90) | Orbital-pair keys for cached Breit Bessel arrays |
| `storage_bytes.f90` | `suite_storage` | [`array_storage_size`](../../src/lib/grasp_shared/array_storage_size.f90) | Reports array storage in bytes |
| `suite_parameters_M.f90` | `suite_parameters` | [`cal_constants`](../../src/lib/grasp_shared/cal_constants.f90) | Shared array capacities, orbital packing radix, and CSF encoding constants |
| `syma_C.f90` | `syma_c` | [`level_symmetry_metadata`](../../src/lib/grasp_shared/level_symmetry_metadata.f90) | Total angular momentum and parity of selected levels |
| `tatb_C.f90` | `tatb_c` | [`radial_quadrature_workspace`](../../src/lib/grasp_shared/radial_quadrature_workspace.f90) | Radial integration work arrays and their active extent |
| `terms_C.f90` | `terms_c` | [`subshell_term_table`](../../src/lib/grasp_shared/subshell_term_table.f90) | Allowed jj-coupled terms indexed by subshell and occupation |
| `titl_C.f90` | `titl_c` | [`calculation_header`](../../src/lib/grasp_shared/calculation_header.f90) | Calculation title and timestamp strings |
| `trk_C.f90` | `trk_c` | [`jj_tensor_coupling_tables`](../../src/lib/grasp_shared/jj_tensor_coupling_tables.f90) | Tensor coupling lookup tables for jj angular algebra |
| `vinlst_C.f90` | `vinlst_c` | [`vinti_integral_table`](../../src/lib/grasp_shared/vinti_integral_table.f90) | Indices and values of Vinti integrals for the specific mass shift |
| `vpilst_C.f90` | `vpilst_c` | [`vacuum_polarization_integral_table`](../../src/lib/grasp_shared/vacuum_polarization_integral_table.f90) | Indices and values of vacuum-polarization radial integrals |
| `wave_C.f90` | `wave_c` | [`radial_wavefunctions`](../../src/lib/grasp_shared/radial_wavefunctions.f90) | Dirac orbital components and radial extent for active and paired bases |
| `wchblk_C.f90` | `wchblk_c` | [`active_symmetry_block`](../../src/lib/grasp_shared/active_symmetry_block.f90) | Index of the symmetry block currently being processed |
| `wfac_C.f90` | `wfac_c` | [`breit_frequency_scale`](../../src/lib/grasp_shared/breit_frequency_scale.f90) | Scale factor for the Breit interaction frequency |
| `whfrom_C.f90` | `whfrom_c` | [`initial_orbital_source`](../../src/lib/grasp_shared/initial_orbital_source.f90) | Source of each initial radial orbital |
