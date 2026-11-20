if(GRASP_BUILD_SERIAL_APPS)
    add_executable(gj
        compute_gj_only_interface.f90
        write_lande_g_factors_interface.f90
        ../rhfs/read_hyperfine_options_interface.f90
        ../rhfs/read_mixing_blocks_interface.f90
        ../rhfs/hyperfine_one_electron_element_interface.f90
        ../rhfs/hyperfine_radial_integral_interface.f90
        ../rhfs/configure_debug_output_interface.f90
        compute_gj_only.f90
        write_lande_g_factors.f90
        ../rhfs/read_hyperfine_options.f90
        ../rhfs/read_mixing_blocks.f90
        ../rhfs/hyperfine_one_electron_element.f90
        ../rhfs/hyperfine_radial_integral.f90
        ../rhfs/configure_debug_output.f90
        gj.f90
    )
    target_link_libraries_Fortran(gj PUBLIC spin_angular breit_angular_support core grasp_shared)
    install(TARGETS gj DESTINATION bin/)
endif()

if(TARGET grasp_mpi)
    add_executable(gj_mpi
        compute_gj_only_mpi_interface.f90
        write_lande_g_factors_interface.f90
        ../rhfs/read_hyperfine_options_interface.f90
        ../rhfs/read_mixing_blocks_interface.f90
        ../rhfs/hyperfine_one_electron_element_interface.f90
        ../rhfs/hyperfine_radial_integral_interface.f90
        ../rhfs/configure_debug_output_interface.f90
        compute_gj_only_mpi.f90
        write_lande_g_factors.f90
        ../rhfs/read_hyperfine_options.f90
        ../rhfs/read_mixing_blocks.f90
        ../rhfs/hyperfine_one_electron_element.f90
        ../rhfs/hyperfine_radial_integral.f90
        ../rhfs/configure_debug_output.f90
        gj_mpi.f90
    )
    target_link_libraries_Fortran(gj_mpi PUBLIC grasp_mpi spin_angular breit_angular_support core grasp_shared)
    target_include_directories(gj_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(gj_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(gj_mpi PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
    )
    install(TARGETS gj_mpi DESTINATION bin/)
endif()
