set(RDENSITY_SOURCES
    read_mixing_blocks.f90
    read_radial_density_input.f90
    configure_debug_output.f90
    natural_orbitals.f90 density_calculation.f90 rdensity.f90
)

if(GRASP_BUILD_SERIAL_APPS)
    add_executable(rdensity density_parallel_serial.f90 ${RDENSITY_SOURCES})
    setup_fortran_modules(rdensity)
    target_link_libraries_Fortran(rdensity PUBLIC spin_angular core grasp_shared)
    install(TARGETS rdensity DESTINATION bin/)
endif()

if(TARGET grasp_mpi)
    add_executable(rdensity_mpi density_parallel_mpi.f90 ${RDENSITY_SOURCES})
    setup_fortran_modules(rdensity_mpi)
    target_link_libraries_Fortran(rdensity_mpi PUBLIC grasp_mpi spin_angular core grasp_shared)
    target_include_directories(rdensity_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(rdensity_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(rdensity_mpi PROPERTIES
        COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
        LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    install(TARGETS rdensity_mpi DESTINATION bin/)
endif()
