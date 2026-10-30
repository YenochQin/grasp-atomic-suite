# Source list follows the GRASPG Makefile (e157d9f).
add_executable(rci_csfg_mpi
    setham_gg.f90
    setham_gg_I.f90
    genmat.f90
    genmat_I.f90
    getcid.f90
    getcid_I.f90
    lodres.f90
    lodres_I.f90
    maneigmpi.f90
    maneigmpi_I.f90
    matrix.f90
    matrix_I.f90
    setres.f90
    setres_I.f90
    rci90mpi.f90
)
setup_fortran_modules(rci_csfg_mpi)
target_include_directories(rci_csfg_mpi PRIVATE "${Fortran_MODULE_DIRECTORY_root}/rci_csfg_mpi")
target_link_libraries_Fortran(rci_csfg_mpi PRIVATE me_csfg dvd90mpi rang90_csfg 9290_csfg mod_csfg mpiu90 9290 mod)
target_link_libraries(rci_csfg_mpi PRIVATE MPI::MPI_Fortran)
target_link_libraries(rci_csfg_mpi PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS rci_csfg_mpi DESTINATION bin/)
