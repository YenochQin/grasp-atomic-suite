# Source list follows the GRASPG Makefile (e157d9f).
add_library(dvd90mpi STATIC
    ../libdvd90/adds_I.f90
    dvdrvr_I.f90
    ../libdvd90/dvdson_I.f90
    ../libdvd90/initdvd_I.f90
    ../libdvd90/mgs_nrm_I.f90
    ../libdvd90/multbc_I.f90
    newvec_I.f90
    ../libdvd90/ovflow_I.f90
    tstsel_I.f90
    dvd_mpi_I.f90
    dvdson.f90
    ../libdvd90/gdvd.f90
    ../libdvd90/gdvd_I.f90
    ../libdvd90/iniest.f90
    ../libdvd90/iniest_I.f90
)
setup_fortran_modules(dvd90mpi)
target_include_directories(dvd90mpi PRIVATE "${Fortran_MODULE_DIRECTORY_root}/dvd90mpi")
target_link_libraries_Fortran(dvd90mpi PRIVATE mpiu90 9290 mod)
target_link_libraries(dvd90mpi PRIVATE MPI::MPI_Fortran)
target_link_libraries(dvd90mpi PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
install(TARGETS dvd90mpi DESTINATION lib/)
