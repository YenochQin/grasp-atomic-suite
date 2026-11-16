add_library(mpiu90 STATIC
    input_io_mpi.f90
    npol_init_mpi.f90
    mpi_C.f90
    mpi_runtime.f90
    mpiu.f90
    cpath.f90
    cslhmpi.f90
    iniestmpi.f90
    lodcslmpi.f90
    lodrwfmpi.f90
    setisompi.f90
    setrwfmpi.f90
    spicmvmpi.f90
    sys_chdir.f90
    sys_getwd.f90
    sys_mkdir.f90
)
setup_fortran_modules(mpiu90)
target_link_libraries_Fortran(mpiu90 PRIVATE mod 9290)
target_include_directories(mpiu90 PRIVATE ${MPI_Fortran_INCLUDE_PATH})
target_link_libraries(mpiu90 PRIVATE ${MPI_Fortran_LIBRARIES})
set_target_properties(mpiu90 PROPERTIES
  COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
install(TARGETS mpiu90 DESTINATION lib/)
# Check NP production routines as well as their test callers in Debug builds.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(npol_init_mpi.f90 PROPERTIES
        COMPILE_OPTIONS "$<$<CONFIG:Debug>:-fcheck=all;-fbacktrace>")
endif()
