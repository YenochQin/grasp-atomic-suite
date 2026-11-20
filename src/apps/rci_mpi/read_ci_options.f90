module read_ci_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_options
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_ci_options (isotope_file_path, radial_file_path, block_labels)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit
        !                                                                      *
        !   Interactively determines the data governing the CI problem.        *
        !   iccut is replaced by an array iccutblk(1:nblock)
        !                                                                      *
        !   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
        !               [RCI92]: SETISO, SETRWF.                               *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
        !   Block version by Xinghong He          Last revision: 15 Jun 1998   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use rci_csf_block_reader, only: read_csf_zero_order_counts
        use mpi_csf_collectives, only: require_collective_success
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
        use standard_io_units, only: stderr_unit
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use nuclear_potential_state,          only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
        use orbital_csf_basis, only: csf_count, orbital_count, orbital_energy_parameters, orbital_origin_exponents
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use breit_frequency_scale, only: transverse_photon_frequency_scale
        use ci_restart_units, only: restart_matrix_unit
        use qed_orbital_cutoffs, only: self_energy_cutoff_flag, self_energy_max_principal_n
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, mpi_logical, mpi_rank
        use mpi_isotope_input, only: read_mpi_isotope_data
        use mpi_radial_wavefunction_file, only: read_mpi_radial_wavefunction_file
        use yes_no_input, only: read_yes_no
        use nuclear_coulomb_potential, only: calculate_nuclear_coulomb_potential
        use radial_grid_setup, only: initialize_radial_grid
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        implicit none (type, external)

        ! Arguments
        character(len=*) :: isotope_file_path, radial_file_path
        character(len=8) ::  block_labels(*)

        ! Local declarations
        logical :: yes
        integer :: i, input_status, np10, orbital_index
        !-----------------------------------------------------------------------
        !
        ! Open, check, load data from, and close the  .iso  file
        ! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole
        ! common block.
        !
        print *, "Calling SETISO ..."
        call read_mpi_isotope_data (isotope_file_path)
        !
        ! The speed of light, if non-default then spread from node-0
        ! Quantities to be obtained: C
        !
        if (nondefault_options_flag /= 0) then
            !cjb node0 input
            if (mpi_rank == 0) then
                write (stderr_unit,*) "Revise the physical speed of light (",vacuum_speed_of_light_au,  &
                    " in a.u.) ?"
                yes = read_yes_no ()
                if (yes) then
                    write (stderr_unit,*) "Enter the revised value:"
                    read *,speed_of_light_au
                else
                    speed_of_light_au = vacuum_speed_of_light_au
                end if
            end if
            call mpi_bcast (speed_of_light_au, 1, mpi_double_precision, 0, &
                mpi_comm_world, mpi_error_code)
        else
            speed_of_light_au = vacuum_speed_of_light_au
        end if
        !
        ! Select the CSF partition independently of the default physical settings.
        input_status = 0
        if (mpi_rank == 0) then
            write (stderr_unit,"(A)") "Use zero-first partition from descriptor header? (y/n)"
            use_first_order_csf_perturbation = read_yes_no ()
            write (calculation_log_unit,"(A)") merge("y", "n", use_first_order_csf_perturbation) &
                & //"            ! Use zero-first partition from descriptor header"
            if (use_first_order_csf_perturbation) call read_csf_zero_order_counts(csf_cutoff_by_block, input_status)
        end if
        call require_collective_success(input_status)
        call mpi_bcast (use_first_order_csf_perturbation, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        if (use_first_order_csf_perturbation) then
            call mpi_bcast (csf_cutoff_by_block, symmetry_block_count, mpi_integer, 0, mpi_comm_world, mpi_error_code)
        else
            csf_cutoff_by_block = csf_counts_by_block(1:symmetry_block_count)
        end if

        !*****************************************************************
        !
        ! Pre-run ?
        !
        !     IF (IPRERUN .EQ. 0) THEN

        !        WRITE (istde,*) ' Prerun with limited interaction?'
        !        YES = GETYN ()

        !        IF (YES) THEN
        !           IPRERUN = 1
        !           LTRANS = .FALSE.
        !           LVP = .FALSE.
        !           LNMS = .FALSE.
        !           LSMS = .FALSE.
        !           LSE = .FALSE.

        !           WRITE (istde,*)  ' Give CSL cut'
        !           READ *, NCSFPRE
        !           WRITE (istde,*)  ' Give coefficient cut for H_0'
        !           READ *, COEFFCUT1
        !           WRITE (istde,*) ' Give coefficient cut for the transvers'
        !    &,                  ' interaction'
        !        ENDIF
        !     ENDIF
        !*****************************************************************
        !
        ! Include transverse ?
        ! Quantities to be obtained: LTRANS, WFACT
        !
        if (mpi_rank == 0) then
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
        end if ! myid .EQ. 0
        call mpi_bcast (include_transverse_interaction, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (transverse_photon_frequency_scale, 1, mpi_double_precision, 0,            &
            mpi_comm_world, mpi_error_code)
        !
        ! Other interactions ? One logical for each case. Done altogether
        !
        if (mpi_rank == 0) then
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
                write (stderr_unit,*)                                           &
                    "Largest n quantum number for including self-energy for orbital"
                write (stderr_unit,*) "n should be less or equal 8"
                read *, self_energy_max_principal_n
            else
                self_energy_cutoff_flag = 0
            end if

            if (include_transverse_interaction) then
                write(calculation_log_unit,"(a)")"y            ! Contribution of H Transverse?"
            else
                write(calculation_log_unit,"(a)")"n            ! Contribution of H Transverse?"
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
        end if ! myid .EQ. 0

        call mpi_bcast (include_vacuum_polarization,  1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (include_normal_mass_shift, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (include_specific_mass_shift, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        call mpi_bcast (estimate_self_energy,  1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
        !
        ! Parameters controlling the radial grid
        !
        !cff-cjb grid
        !cff use default grid
        !
        ! 99  IF (NPARM .EQ. 0) THEN
        !        RNT = EXP (-65.D0/16.D0) / Z
        !        H = 0.5D0**4
        !        N = MIN (220,NNNP)
        !     ELSE
        !CFF     .. should be Z-dependent
        radial_grid_scale = 2.0e-6_real64/nuclear_charge
        transformed_grid_step = 5.0e-2_real64
        radial_point_count = max_radial_points
        !     ENDIF
        asymptotic_grid_step = 0.0_real64
        !     IF ( NDEF.NE.0) THEN
        !        IF (myid .EQ. 0) THEN
        !           WRITE (istde,*) 'The default radial grid parameters',      &
        !                          ' for this case are:'
        !           WRITE (istde,*) ' RNT = ',RNT,';'
        !           WRITE (istde,*) ' H = ',H,';'
        !           WRITE (istde,*) ' HP = ',HP,';'
        !           WRITE (istde,*) ' N = ',N,';'
        !           WRITE (istde,*) ' revise these values?'
        !           YES = GETYN ()
        !        ENDIF
        !
        !         ...To prevent subsequent BCAST when YES is false
        !        CALL MPI_Bcast (YES, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        !        IF (YES) THEN
        !           IF (myid .EQ. 0) THEN
        !              WRITE (istde,*) 'Enter RNT:'
        !              READ *, RNT
        !              WRITE (istde,*) 'Enter H:'
        !              READ *, H
        !              WRITE (istde,*) 'Enter HP:'
        !              READ *, HP
        !              WRITE (istde,*) 'Enter N:'
        !              READ *, N
        !           ENDIF
        !           CALL MPI_Bcast (RNT, 1, MPI_DOUBLE_PRECISION, 0,    &
        !                               MPI_COMM_WORLD, ierr)
        !           CALL MPI_Bcast (H,   1, MPI_DOUBLE_PRECISION, 0,    &
        !                               MPI_COMM_WORLD, ierr)
        !           CALL MPI_Bcast (HP,  1, MPI_DOUBLE_PRECISION, 0,    &
        !                               MPI_COMM_WORLD, ierr)
        !           CALL MPI_Bcast (N,   1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        !        ENDIF
        !     ENDIF ! NDEF.NE.0
        !
        !   ACCY is an estimate of the accuracy of the numerical procedures
        !
        working_tolerance = transformed_grid_step**6
        !
        !   Set up the coefficients for the numerical procedures
        !
        call initialize_radial_stencil_weights
        !
        !   Generate the radial grid and all associated arrays
        !
        call initialize_radial_grid
        !
        !   Generate $- r \times V_ (r)$
        !
        call calculate_nuclear_coulomb_potential
        !
        !   Load the radial wavefunctions
        !
        call read_mpi_radial_wavefunction_file (radial_file_path)
        !
        !   Write the basic parameters of the model electron cloud to the
        !   .res  file; this is the second record on the file --- the
        !   first is the header (see SUBROUTINE SETRES)
        !
        write (restart_matrix_unit) electron_count, csf_count, orbital_count, symmetry_block_count ! ncf is ncftot, from setcsll
        !
        !   Write the nuclear parameters and $- r \times V_ (r)$
        !   to the  .res  file; these are the third, fourth, and fifth
        !   records on the file
        !
        write (restart_matrix_unit) nuclear_charge, nuclear_mass_in_electron_masses
        write (restart_matrix_unit) nuclear_model_parameter_count,(nuclear_model_parameters(i), i = 1, &
            & nuclear_model_parameter_count)
        write (restart_matrix_unit) radial_point_count, (negative_radius_times_potential(i), i = 1, radial_point_count), &
            & first_grid_index_beyond_fermi_c
        !
        !   Write the physical effects specification to the  .res  file.
        !   iccutblk() is now an array of length nblock.
        !   This is the sixth record on the file
        !
        write (restart_matrix_unit) speed_of_light_au, use_first_order_csf_perturbation, (csf_cutoff_by_block(i), i = 1, &
            symmetry_block_count),           &
            include_transverse_interaction, transverse_photon_frequency_scale, include_vacuum_polarization, &
                include_normal_mass_shift, include_specific_mass_shift
        !
        !   Write the grid data to the  .res  file; this is the seventh
        !   record on the file
        !
        np10 = radial_point_count + 10
        write (restart_matrix_unit) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, (radial_grid(i), i = 1, np10), &
            (radial_grid_derivative(i), i = 1, np10), (radial_derivative_over_radius(i), i = 1, np10) ! (restart_matrix_unit)
        !
        !   Write out the interpolated radial wavefunctions; there are
        !   2*NW such records; thus the total number of records written
        !   at this point is 7+2*NW
        !
        do orbital_index = 1, orbital_count
            write (restart_matrix_unit) orbital_energy_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
                & orbital_origin_parameters(orbital_index), &
                & orbital_last_grid_index(orbital_index)
            write (restart_matrix_unit) (radial_large_components(i,orbital_index), i = 1, orbital_last_grid_index(orbital_index)), &
                (radial_small_components(i,orbital_index) &
                & , i = 1, orbital_last_grid_index(orbital_index))
        end do
        !
        return
    end subroutine read_ci_options
end module read_ci_options_interface
