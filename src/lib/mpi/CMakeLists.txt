add_library(grasp_mpi STATIC
    mpi_csf_block_input.f90
    mpi_csf_collectives.f90
    mpi_csf_basis_input.f90
    mpi_rmcdhf_csf_input.f90
    mpi_io_errors.f90
    mpi_nuclear_polarization.f90
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
target_include_directories(grasp_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(grasp_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
set_target_properties(grasp_mpi PROPERTIES
  COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
install(TARGETS grasp_mpi DESTINATION lib/)
# Check NP production routines as well as their test callers in Debug builds.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(mpi_nuclear_polarization.f90 PROPERTIES
        COMPILE_OPTIONS "$<$<CONFIG:Debug>:-fcheck=all;-fbacktrace>")
endif()
