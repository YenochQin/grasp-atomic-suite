if(GRASP_BUILD_SERIAL_APPS)
    add_executable(hyperfine
        write_hyperfine_energies.f90
        read_hyperfine_options.f90
        read_mixing_blocks.f90
        hyperfine_one_electron_element.f90
        hyperfine_radial_integral.f90
        hyperfine_nuclear_factors.f90
        configure_debug_output.f90
        open_summary_file.f90
        write_input_summary.f90
        calculate_hyperfine_matrix_elements.f90
        hyperfine_main.f90
)
    target_link_libraries_Fortran(hyperfine PUBLIC spin_angular breit_angular_support core grasp_shared)
    install(TARGETS hyperfine DESTINATION bin/)
endif()

if(TARGET grasp_mpi)
    add_executable(hyperfine_mpi
        write_hyperfine_energies.f90
        read_hyperfine_options.f90
        read_mixing_blocks.f90
        hyperfine_one_electron_element.f90
        hyperfine_radial_integral.f90
        hyperfine_nuclear_factors.f90
        configure_debug_output.f90
        open_summary_file.f90
        write_input_summary.f90
        calculate_hyperfine_matrix_elements_mpi.f90
        hyperfine_mpi_main.f90
)
    target_link_libraries_Fortran(hyperfine_mpi PUBLIC grasp_mpi spin_angular breit_angular_support core grasp_shared)
    install(TARGETS hyperfine_mpi DESTINATION bin/)
endif()
