add_executable(rwfnestimate
    rwfnestimate_input.f90
    integrate_radial_tail.f90
    integrate_radial_tail_interface.f90
    integrate_homogeneous_radial_segment.f90
    integrate_homogeneous_radial_segment_interface.f90
    solve_homogeneous_radial_equation.f90
    solve_homogeneous_radial_equation_interface.f90
    estimate_hydrogenic_orbitals.f90
    estimate_hydrogenic_orbitals_interface.f90
    estimate_from_saved_orbitals.f90
    estimate_from_saved_orbitals_interface.f90
    estimate_thomas_fermi_orbitals.f90
    estimate_thomas_fermi_orbitals_interface.f90
    build_thomas_fermi_potential.f90
    build_thomas_fermi_potential_interface.f90
    print_uninitialized_orbitals.f90
    print_uninitialized_orbitals_interface.f90
    set_hydrogenic_screening.f90
    set_hydrogenic_screening_interface.f90
    configure_debug_output.f90
    configure_debug_output_interface.f90
    open_summary_file.f90
    open_summary_file_interface.f90
    write_input_summary.f90
    write_input_summary_interface.f90
    write_initial_wavefunction_summary.f90
    write_initial_wavefunction_summary_interface.f90
    write_initial_radial_wavefunctions.f90
    write_initial_radial_wavefunctions_interface.f90
    generate_initial_radial_wavefunctions.f90
    generate_initial_radial_wavefunctions_interface.f90
    read_calculation_options.f90
    read_calculation_options_interface.f90
    rwfnestimate_main.f90
)
target_link_libraries_Fortran(rwfnestimate PUBLIC spin_angular core grasp_shared)

# Keep the header parser and discard the unused Parquet code in the Rust archive.
include(CheckLinkerFlag)
if(APPLE)
    set(rwfnestimate_dead_strip "LINKER:-dead_strip")
else()
    set(rwfnestimate_dead_strip "LINKER:--gc-sections")
endif()
check_linker_flag(Fortran "${rwfnestimate_dead_strip}" GRASP_RWFNESTIMATE_DEAD_STRIP)
if(GRASP_RWFNESTIMATE_DEAD_STRIP)
    target_link_options(rwfnestimate PRIVATE "${rwfnestimate_dead_strip}")
endif()

install(TARGETS rwfnestimate DESTINATION bin/)
