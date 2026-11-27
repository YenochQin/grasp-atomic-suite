if(GRASP_BUILD_SERIAL_APPS)
    add_executable(gj

        compute_gj_only.f90
        write_lande_g_factors.f90
        ../hyperfine/read_hyperfine_options.f90
        ../hyperfine/read_mixing_blocks.f90
        ../hyperfine/hyperfine_one_electron_element.f90
        ../hyperfine/hyperfine_radial_integral.f90
        ../hyperfine/configure_debug_output.f90
        gj.f90
)
    target_link_libraries_Fortran(gj PUBLIC spin_angular breit_angular_support core grasp_shared)
    install(TARGETS gj DESTINATION bin/)
endif()

if(TARGET grasp_mpi)
    add_executable(gj_mpi

        compute_gj_only_mpi.f90
        write_lande_g_factors.f90
        ../hyperfine/read_hyperfine_options.f90
        ../hyperfine/read_mixing_blocks.f90
        ../hyperfine/hyperfine_one_electron_element.f90
        ../hyperfine/hyperfine_radial_integral.f90
        ../hyperfine/configure_debug_output.f90
        gj_mpi.f90
)
    target_link_libraries_Fortran(gj_mpi PUBLIC grasp_mpi spin_angular breit_angular_support core grasp_shared)
    install(TARGETS gj_mpi DESTINATION bin/)
endif()
