add_library(grasp_mpi STATIC
    mpi_csf_block_input.f90
    mpi_csf_collectives.f90
    mpi_csf_basis_input.f90
    mpi_rmcdhf_csf_input.f90
    mpi_io_errors.f90
    mpi_element_counts.f90
    mpi_indexed_mcp.f90
    mpi_mcp_calculation_input.f90
    mpi_state.f90
    mpi_runtime.f90
    mpi_run_logging.f90
    mpi_utilities.f90
    mpi_work_paths.f90
    mpi_calculation_state.f90
    mpi_initial_eigenvectors.f90
    mpi_radial_wavefunction_reader.f90
    mpi_isotope_input.f90
    mpi_radial_wavefunction_file.f90
    mpi_sparse_matrix_product.f90
    mpi_change_directory.f90
    mpi_current_directory.f90
    mpi_make_directory.f90
)
setup_fortran_modules(grasp_mpi)
target_link_libraries_Fortran(grasp_mpi PRIVATE grasp_shared core)
target_link_libraries(grasp_mpi PUBLIC MPI::MPI_Fortran)
install(TARGETS grasp_mpi DESTINATION lib/)
