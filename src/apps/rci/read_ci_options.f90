! Set up the nuclear and radial data and choose the interactions included in the CI Hamiltonian.
module read_ci_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_options
contains

    !
    ! Written by Farid A. Parpia. Last revision: 15 Dec 1992
    ! Block version by Xinghong He. Last revision: 15 Jun 1998
    subroutine read_ci_options (isotope_file_path, radial_file_path, block_labels)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation, include_normal_mass_shift, &
            estimate_self_energy, include_specific_mass_shift, include_transverse_interaction, include_vacuum_polarization
        use atomic_calculation_state,           only: vacuum_speed_of_light_au, speed_of_light_au, nuclear_charge, &
            working_tolerance, electron_count, nuclear_mass_in_electron_masses
        use calculation_input_defaults, only: nondefault_options_flag
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_derivative, radial_derivative_over_radius, &
            radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
        use csf_block_cutoffs, only: csf_cutoff_by_block
        use standard_io_units, only: stderr_unit, stdin_unit
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state,          only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_energy_parameters, orbital_origin_exponents
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use breit_frequency_scale, only: transverse_photon_frequency_scale
        use ci_restart_units, only: restart_matrix_unit
        use qed_orbital_cutoffs, only: self_energy_cutoff_flag, self_energy_max_principal_n
        use yes_no_input, only: read_yes_no
        use nuclear_coulomb_potential, only: calculate_nuclear_coulomb_potential
        use radial_grid_setup, only: initialize_radial_grid
        use isotope_data_input, only: read_isotope_data
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_wavefunction_input, only: read_radial_wavefunction_file
        implicit none (type, external)

        ! Arguments
        character(len=*) :: isotope_file_path, radial_file_path
        character(len=8) ::  block_labels(*)

        ! Local declarations
        logical :: yes
        integer :: symmetry_block_index, jblock, ntmp, np10, orbital_index
        ! Load the nuclear and isotope parameters for the CI Hamiltonian.
        call read_isotope_data (isotope_file_path)
        ! Allow the default speed of light to be overridden for this calculation.
        if (nondefault_options_flag /= 0) then
            write (stderr_unit,*) "Revise the physical speed of light (",vacuum_speed_of_light_au,   &
                " in a.u.) ?"
            yes = read_yes_no ()
            if (yes) then
                write (stderr_unit,*) "Enter the revised value:"
                read *,speed_of_light_au
            else
                speed_of_light_au = vacuum_speed_of_light_au
            end if
        else
            speed_of_light_au = vacuum_speed_of_light_au
        end if
        ! Select any first-order partition of the CSF basis.
        if (nondefault_options_flag /= 0) then
            write (stderr_unit,*) "Treat contributions of some CSFs",            &
                " as first-order perturbations?"
            use_first_order_csf_perturbation = read_yes_no ()
        else
            use_first_order_csf_perturbation = .false.
        end if

        ! Set the cutoff for first-order CSF interactions in each symmetry block.

        if (.not. use_first_order_csf_perturbation) then
            ! Use the entire block when first-order restrictions are disabled.
            do symmetry_block_index = 1, symmetry_block_count
                csf_cutoff_by_block(symmetry_block_index) = csf_counts_by_block(symmetry_block_index)
            end do
        else

            ! Read the first-order cutoff for each symmetry block.
            write (stderr_unit,*) "There are ", symmetry_block_count, "blocks. They are:"
            write (stderr_unit,*) "  block     J Parity     No of CSFs"
            do symmetry_block_index = 1, symmetry_block_count
                write (stderr_unit,*) symmetry_block_index, trim(block_labels(symmetry_block_index)), &
                    & csf_counts_by_block(symmetry_block_index)
            end do

            write (stderr_unit,*)
            write (stderr_unit,*) "Enter iccut for each block"
            do jblock = 1, symmetry_block_count
                write (stderr_unit,*) "Block ", jblock, "   ncf = ",           &
                    csf_counts_by_block(jblock)                            &
                    , " id = ", trim(block_labels(jblock))
                do
                    read (stdin_unit,*) ntmp
                    if (ntmp >= 0 .and. ntmp <= csf_counts_by_block(jblock)) exit
                    write (stderr_unit,*) "ICCUT out of range, re-enter:"
                end do
                csf_cutoff_by_block(jblock) = ntmp
                write(calculation_log_unit,*) ntmp,"! ICCUT for block",jblock
            end do
        end if

        ! Select the Breit interaction and optional frequency scaling.
        write (stderr_unit,*) "Include contribution of H (Transverse)?"
        include_transverse_interaction = read_yes_no ()
        write (stderr_unit,*) "Modify all transverse photon frequencies?"
        yes = read_yes_no ()
        if (yes) then
            write (stderr_unit,*) "Enter the scale factor:"
            read *, transverse_photon_frequency_scale
        else
            transverse_photon_frequency_scale = 1.0_real64
        end if
        ! Select the remaining mass-shift and QED corrections.
        write (stderr_unit,*) "Include H (Vacuum Polarisation)?"
        include_vacuum_polarization = read_yes_no ()

        write (stderr_unit,*) "Include H (Normal Mass Shift)?"
        include_normal_mass_shift = read_yes_no ()

        write (stderr_unit,*) "Include H (Specific Mass Shift)?"
        include_specific_mass_shift = read_yes_no ()

        write (stderr_unit,*) "Estimate self-energy?"
        estimate_self_energy = read_yes_no ()
        if (estimate_self_energy.eqv..true.) then
            self_energy_cutoff_flag = 1
            write (stderr_unit,*)                                             &
                "Largest n quantum number for including self-energy for orbital"
            write (stderr_unit,*) "n should be less or equal 8"
            read *, self_energy_max_principal_n
        else
            self_energy_cutoff_flag = 0
        end if

        if (include_transverse_interaction) then
            write(calculation_log_unit,"(a)") "y            ! Contribution of H Transverse?"
        else
            write(calculation_log_unit,"(a)") "n            ! Contribution of H Transverse?"
        end if
        if (yes) then
            write(calculation_log_unit,"(a)") "y            ! Modify photon frequencies?"
            write(calculation_log_unit,*) transverse_photon_frequency_scale,"! Scale factor"
        else
            write(calculation_log_unit,"(a)") "n            ! Modify photon frequencies?"
        end if
        if (include_vacuum_polarization) then
            write(calculation_log_unit,"(a)") "y            ! Vacuum polarization?"
        else
            write(calculation_log_unit,"(a)") "n            ! Vacuum polarization?"
        end if
        if (include_normal_mass_shift) then
            write(calculation_log_unit,"(a)") "y            ! Normal mass shift?"
        else
            write(calculation_log_unit,"(a)") "n            ! Normal mass shift?"
        end if
        if (include_specific_mass_shift) then
            write(calculation_log_unit,"(a)") "y            ! Specific mass shift?"
        else
            write(calculation_log_unit,"(a)") "n            ! Specific mass shift?"
        end if
        if (estimate_self_energy) then
            write(calculation_log_unit,"(a)") "y            ! Self energy?"
            write(calculation_log_unit,*) self_energy_max_principal_n, "! Max n for including self energy"
        else
            write(calculation_log_unit,"(a)") "n            ! Self energy?"
        end if
        !  Parameters controlling the radial grid
        !
        ! cjb grid
        ! cff use default grid
        !
        !         RNT = EXP (-65.D0/16.D0) / Z
        !         H = 0.5D0**4
        !         N = MIN (220,NNNP)
        ! The screening estimate could be refined to depend on nuclear charge Z.
        radial_grid_scale = 2.0e-6_real64/nuclear_charge
        transformed_grid_step = 5.0e-2_real64
        radial_point_count = max_radial_points
        asymptotic_grid_step = 0.0_real64
        ! ACCY is an estimate of the accuracy of the numerical procedures
        working_tolerance = transformed_grid_step**6
        ! Initialize the radial integration and derivative stencils.
        call initialize_radial_stencil_weights
        ! Build the radial grid and its derivative arrays.
        call initialize_radial_grid
        ! Generate $- r \times V_ (r)$
        call calculate_nuclear_coulomb_potential
        ! Load the radial wavefunctions
        call read_radial_wavefunction_file (radial_file_path)
        ! Write the basic parameters of the model electron cloud to the
        ! .res  file; this is the second record on the file --- the
        ! first is the header (see SUBROUTINE SETRES)
        write (restart_matrix_unit) electron_count, csf_count, orbital_count, symmetry_block_count ! ncf is ncftot, from setcsll
        ! Write the nuclear parameters and $- r \times V_ (r)$
        ! to the  .res  file; these are the third, fourth, and fifth
        ! records on the file
        write (restart_matrix_unit) nuclear_charge, nuclear_mass_in_electron_masses
        write (restart_matrix_unit) nuclear_model_parameter_count,(nuclear_model_parameters(symmetry_block_index), &
            & symmetry_block_index = 1, nuclear_model_parameter_count)
        write (restart_matrix_unit) radial_point_count, (negative_radius_times_potential(symmetry_block_index), &
            & symmetry_block_index = 1, radial_point_count), &
            & first_grid_index_beyond_fermi_c
        ! Write the physical effects specification to the  .res  file.
        ! iccutblk() is now an array of length nblock.
        ! This is the sixth record on the file
        write (restart_matrix_unit) speed_of_light_au, use_first_order_csf_perturbation, &
            (csf_cutoff_by_block(symmetry_block_index), &
            symmetry_block_index = 1, &
            & symmetry_block_count),            &
            include_transverse_interaction, transverse_photon_frequency_scale, include_vacuum_polarization, &
                include_normal_mass_shift, include_specific_mass_shift
        ! Write the grid data to the  .res  file; this is the seventh
        ! record on the file
        np10 = radial_point_count + 10
        write (restart_matrix_unit) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, &
            & (radial_grid(symmetry_block_index), symmetry_block_index = 1, np10), &
            (radial_grid_derivative(symmetry_block_index), symmetry_block_index = 1, np10), &
            & (radial_derivative_over_radius(symmetry_block_index), symmetry_block_index = 1, np10) ! (restart_matrix_unit)
        ! Write out the interpolated radial wavefunctions; there are
        ! 2*NW such records; thus the total number of records written
        ! at this point is 7+2*NW
        do orbital_index = 1, orbital_count
            write (restart_matrix_unit) orbital_energy_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
                & orbital_origin_parameters(orbital_index), &
                & orbital_last_grid_index(orbital_index)
            write (restart_matrix_unit) (radial_large_components(symmetry_block_index,orbital_index), symmetry_block_index = 1, &
                & orbital_last_grid_index(orbital_index)), &
                (radial_small_components(symmetry_block_index,orbital_index) &
                & , symmetry_block_index = 1, orbital_last_grid_index(orbital_index))
        end do
        return
    end subroutine read_ci_options
end module read_ci_options_interface
