if(GRASP_BUILD_SERIAL_APPS)
    add_executable(rhfs
        write_hyperfine_energies.f90
        read_hyperfine_options.f90
        read_mixing_blocks.f90
        hyperfine_one_electron_element.f90
        hyperfine_radial_integral.f90
        configure_debug_output.f90
        open_summary_file.f90
        write_input_summary.f90
        calculate_hyperfine_matrix_elements.f90
        rhfs_main.f90
)
    target_link_libraries_Fortran(rhfs PUBLIC spin_angular breit_angular_support core grasp_shared)
    install(TARGETS rhfs DESTINATION bin/)
endif()

if(TARGET grasp_mpi)
    add_executable(rhfs_mpi
        write_hyperfine_energies.f90
        read_hyperfine_options.f90
        read_mixing_blocks.f90
        hyperfine_one_electron_element.f90
        hyperfine_radial_integral.f90
        configure_debug_output.f90
        open_summary_file.f90
        write_input_summary.f90
        calculate_hyperfine_matrix_elements_mpi.f90
        rhfs_mpi_main.f90
)
    target_link_libraries_Fortran(rhfs_mpi PUBLIC grasp_mpi spin_angular breit_angular_support core grasp_shared)
    target_include_directories(rhfs_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(rhfs_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(rhfs_mpi PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
    install(TARGETS rhfs_mpi DESTINATION bin/)
endif()
