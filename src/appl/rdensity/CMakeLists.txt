set(RDENSITY_SOURCES
    getmixblock.f90 getmixblock_I.f90
    getsmd.f90 getsmd_I.f90
    setdbg.f90 setdbg_I.f90
    natural_orbitals.f90 natorbnew.f90 rdensity.f90
)

if(GRASP_BUILD_SERIAL_APPS)
    add_executable(rdensity density_parallel_serial.f90 ${RDENSITY_SOURCES})
    setup_fortran_modules(rdensity)
    target_link_libraries_Fortran(rdensity PUBLIC rang90 9290 mod)
    install(TARGETS rdensity DESTINATION bin/)
endif()

if(TARGET mpiu90)
    add_executable(rdensity_mpi density_parallel_mpi.f90 ${RDENSITY_SOURCES})
    setup_fortran_modules(rdensity_mpi)
    target_link_libraries_Fortran(rdensity_mpi PUBLIC mpiu90 rang90 9290 mod)
    target_include_directories(rdensity_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(rdensity_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(rdensity_mpi PROPERTIES
        COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
        LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    install(TARGETS rdensity_mpi DESTINATION bin/)
endif()
