! Restore the CI calculation controls and radial state from the restart file.
module read_ci_restart_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_restart
contains

    ! Loads the data from the  .res  file. A number of checks are made
    ! to ensure correctness and consistency.
    !
    ! Written by Farid A. Parpia. Last revision: 15 Oct 1992
    ! Block version by Xinghong He. Last revision:  1 Jun 1998
    subroutine read_ci_restart
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: calculation_log_unit
        use grasp_parameters,   only: max_radial_points
        use array_memory_management, only: alloc
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation, include_normal_mass_shift, &
            estimate_self_energy, include_specific_mass_shift, include_transverse_interaction, include_vacuum_polarization
        use atomic_calculation_state, only: left_electron_count, electron_count, nuclear_charge, &
            nuclear_mass_in_electron_masses, speed_of_light_au, working_tolerance
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_derivative, radial_derivative_over_radius, &
            radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state, only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_energy_parameters, orbital_origin_exponents
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use breit_frequency_scale, only: transverse_photon_frequency_scale
        use ci_restart_units, only: restart_matrix_unit
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        use csf_block_cutoffs, only: csf_cutoff_by_block
        use standard_io_units, only: stderr_unit
        use yes_no_input, only: read_yes_no
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        implicit none (type, external)

        ! Local declarations
        integer :: ncfres, nwres, nblockres, i, np10, orbital_index

        write (output_unit, *) "Calling read_ci_restart ..."
        ! Read the basic parameters of the electron cloud; check these
        ! against the invariants recorded in the descriptor basis.
        read (restart_matrix_unit) left_electron_count, ncfres, nwres, nblockres

        if (left_electron_count/=electron_count .or. ncfres/=csf_count .or. nwres/=orbital_count .or. &
            nblockres/=symmetry_block_count&
            ) stop "read_ci_restart: electron, CSF, or orbital count mismatch"
        ! Read the nuclear parameters
        read (restart_matrix_unit) nuclear_charge, nuclear_mass_in_electron_masses
        read (restart_matrix_unit) nuclear_model_parameter_count, (nuclear_model_parameters(i),i=1,nuclear_model_parameter_count)
        read (restart_matrix_unit) radial_point_count, (negative_radius_times_potential(i),i=1,radial_point_count), &
            first_grid_index_beyond_fermi_c

        if (radial_point_count > max_radial_points) stop "read_ci_restart: radial point count exceeds max_radial_points"
        ! Read the physical effects specifications
        ! iccutblk() is now an array of length nblock.
        read (restart_matrix_unit) speed_of_light_au, use_first_order_csf_perturbation, (csf_cutoff_by_block(i),i=1, &
            symmetry_block_count), include_transverse_interaction, transverse_photon_frequency_scale, include_vacuum_polarization, &
            include_normal_mass_shift, include_specific_mass_shift
        ! Read the remaining parameters controlling the radial grid and the
        ! grid arrays
        np10 = radial_point_count + 10
        read (restart_matrix_unit) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, (radial_grid(i),i=1,np10), &
            & (radial_grid_derivative(i),i=1,np10), (radial_derivative_over_radius(i),i=1,&
            np10)
        ! ACCY is an estimate of the accuracy of the numerical procedures
        working_tolerance = transformed_grid_step**6
        ! Initialize the radial integration and derivative stencils.
        call initialize_radial_stencil_weights
        ! Allocate storage for the radial wavefunction arrays
        call alloc (radial_large_components, max_radial_points,orbital_count, "radial_large_components", "read_ci_restart")
        call alloc (radial_small_components, max_radial_points,orbital_count, "radial_small_components", "read_ci_restart")
        ! Read the orbital wavefunctions and the associated arrays
        do orbital_index = 1, orbital_count
            read (restart_matrix_unit) orbital_energy_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
                & orbital_origin_parameters(orbital_index), &
                & orbital_last_grid_index(orbital_index)
            read (restart_matrix_unit) (radial_large_components(i,orbital_index),i=1,orbital_last_grid_index(orbital_index)), &
                & (radial_small_components(i, &
                orbital_index),i=1, &
                & orbital_last_grid_index(orbital_index))
        end do
        ! Determine if the self-energy contribution is to be estimated
        write (stderr_unit, *) "Estimate contributions from the self-energy?"
        estimate_self_energy = read_yes_no()
        if (estimate_self_energy) then
            write(calculation_log_unit,"(a)") "y            ! Estimate contributions from the self-energy?"
        else
            write(calculation_log_unit,"(a)") "n            ! Estimate contributions from the self-energy?"
        end if

        return
    end subroutine read_ci_restart
end module read_ci_restart_interface
