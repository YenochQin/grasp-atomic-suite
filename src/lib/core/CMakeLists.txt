add_library(core STATIC
    calculation_state_files.f90
    csf_descriptor_c_bindings.f90
    rmcdhf_csf_input.f90
    csf_full_basis_reader.f90
    csf_basis_merging.f90
    rci_csf_block_reader.f90
    rci_restart_records.f90
    csf_descriptor_reader.f90
    csf_quantum_queries.f90
    posix_services.f90
    nuclear_polarization_input.f90
    nuclear_polarization_potential.f90
    nuclear_polarization_initialization.f90
    angular_coefficient_buffer.f90
    nonnegative_atan2.f90
    calendar_clock.f90
    complex_gamma.f90
    relativistic_three_j.f90
    integer_format.f90
    integer_format_with_context.f90
    half_integer_format.f90
    mcp_coefficient_output.f90
    closed_shell_mcp.f90
    radial_node_count.f90
    reduced_angular_matrix_element.f90
    file_open.f90
    dirac_coulomb_radial_orbital.f90
    real_vector_initialization.f90
    sparse_dot_and_scatter.f90
    radial_wavefunction_derivatives.f90
    radial_ascii_plot.f90
    alternating_exponential_series.f90
    log_factorial_table.f90
    decimal_digit_check.f90
    subshell_label_selection.f90
    yes_no_input.f90
    one_electron_csf_filter.f90
    two_electron_csf_filter.f90
    subshell_occupation_class.f90
    integer_vector_copy.f90
    sparse_eigensolver_initial_basis.f90
    aitken_interpolation.f90
    radial_quadrature.f90
    weighted_radial_overlap.f90
    radial_orbital_interpolation.f90
    level_selection_parser.f90
    angular_momentum_triangle.f90
    csf_byte_packing.f90
    isotope_data_records.f90
    radial_orbital_orthogonalization.f90
    mixing_file_format.f90
    radial_wavefunction_format.f90
    radial_wavefunction_records.f90
    block_level_selection.f90
    subshell_parent_table.f90
    nuclear_coulomb_potential.f90
    file_positioning.f90
    radial_grid_setup.f90
    one_electron_radial_integral.f90
    physical_constants.f90
    csf_basis_loader.f90
    isotope_data_input.f90
    angular_recoupling_tables.f90
    machine_constants.f90
    radial_dirac_coefficients.f90
    radial_stencil_weights.f90
    csf_pair_quantum_numbers.f90
    radial_wavefunction_input.f90
    coulomb_tensor_rank_range.f90
    hartree_yz_functions.f90
    slater_radial_integral.f90
    sparse_symmetric_matvec.f90
    radial_origin_solution.f90
    run_timer.f90
)
setup_fortran_modules(core)
target_link_libraries(core PUBLIC csf_descriptor)
target_link_libraries_Fortran(core PRIVATE grasp_shared)
install(TARGETS core DESTINATION lib/)
# Enable GNU Fortran runtime checks and backtraces for the nuclear-polarization routines in Debug builds.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(nuclear_polarization_input.f90 nuclear_polarization_potential.f90 nuclear_polarization_initialization.f90 PROPERTIES
        COMPILE_OPTIONS "$<$<CONFIG:Debug>:-fcheck=all;-fbacktrace>")
endif()
