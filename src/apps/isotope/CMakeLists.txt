add_executable(fical isotope_shift_transition_main.f90)
target_link_libraries_Fortran(fical PRIVATE core)
install(TARGETS fical DESTINATION bin/)

if(GRASP_BUILD_SERIAL_APPS)
    add_executable(isotope
        isotope_shift_options.f90
        isotope_main.f90
        open_isotope_shift_angular_cache.f90
        read_field_shift_factors.f90
        read_seltzer_field_shift_factors.f90
        calculate_field_shift_factors.f90
        calculate_seltzer_field_shift_factors.f90
        fit_electron_density.f90
        read_mixing_blocks.f90
        prepare_isotope_calculation.f90
        interpolate_electron_density.f90
        calculate_isotope_shifts.f90
        field_shift_radial_integral.f90
        field_shift_radial_integral_vector.f90
        normal_mass_shift_radial_integral.f90
        specific_mass_shift_radial_integral_2.f90
        specific_mass_shift_radial_integral_3.f90
        configure_debug_output.f90
        open_summary_file.f90
        write_input_summary.f90
        read_specific_mass_shift.f90
        calculate_specific_mass_shift.f90
        first_spin_angular_factor.f90
        second_spin_angular_factor.f90
        vinti_radial_integral.f90
        write_level_energies.f90
)
    target_link_libraries_Fortran(isotope PUBLIC grasp_shared core spin_angular breit_angular_support grasp_shared core)
    install(TARGETS isotope DESTINATION bin/)
endif()

if(TARGET grasp_mpi)
    add_executable(isotope_mpi
        isotope_shift_options.f90
        open_isotope_shift_angular_cache.f90
        read_field_shift_factors.f90
        read_field_shift_factors_mpi.f90
        read_seltzer_field_shift_factors.f90
        read_seltzer_field_shift_factors_mpi.f90
        calculate_field_shift_factors.f90
        calculate_field_shift_factors_mpi.f90
        calculate_seltzer_field_shift_factors.f90
        calculate_seltzer_field_shift_factors_mpi.f90
        fit_electron_density.f90
        read_mixing_blocks.f90
        prepare_isotope_calculation.f90
        interpolate_electron_density.f90
        calculate_isotope_shifts_mpi.f90
        field_shift_radial_integral.f90
        field_shift_radial_integral_vector.f90
        normal_mass_shift_radial_integral.f90
        specific_mass_shift_radial_integral_2.f90
        specific_mass_shift_radial_integral_3.f90
        configure_debug_output.f90
        open_summary_file.f90
        calculate_specific_mass_shift.f90
        read_specific_mass_shift.f90
        read_specific_mass_shift_mpi.f90
        calculate_specific_mass_shift_mpi.f90
        first_spin_angular_factor.f90
        second_spin_angular_factor.f90
        vinti_radial_integral.f90
        write_level_energies.f90
        isotope_mpi_main.f90
)
    target_link_libraries_Fortran(isotope_mpi PUBLIC grasp_mpi grasp_shared core spin_angular breit_angular_support grasp_shared core)
    install(TARGETS isotope_mpi DESTINATION bin/)
endif()
