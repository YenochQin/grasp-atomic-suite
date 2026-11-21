module generate_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: generate_initial_radial_wavefunctions
contains

    !***********************************************************************
    !                                                                      *
    subroutine generate_initial_radial_wavefunctions(default_orbitals)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit, summary_output_unit
        !                                                                      *
        !   Controls the computation of the subshell radial wavefunctions.     *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, GETRSL, GETYN.                         *
        !               [ERWF]: FRMHYD, FRMRWF, FRMTFP, PRTREM, SUMMRY,        *
        !                       TFPOT.                                         *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
        !   Updated by Xinghong He                Last revision: 11 Nov 1997   *
        !    for menu-driven and less questions
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc
        use radial_node_count_settings, only: relative_extremum_cutoff
        use atomic_calculation_state, only: nuclear_charge, speed_of_light_au, max_radial_points
        use calculation_input_defaults, only: nondefault_options_flag
        use radial_grid_state, only: radial_point_count
        use initial_orbital_selection, only: orbital_initialized
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_origin_exponents, &
            orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions, only: radial_large_components, radial_small_components
        use initial_orbital_source, only: orbital_source_description
        use standard_io_units, only: stderr_unit, stdin_unit
        use yes_no_input, only: read_yes_no

        ! Procedure interfaces
        use build_thomas_fermi_potential_interface, only: build_thomas_fermi_potential
        use print_uninitialized_orbitals_interface, only: print_uninitialized_orbitals

        ! Module dependencies
        use subshell_label_selection, only: select_subshell_labels

        ! Procedure interfaces
        use estimate_from_saved_orbitals_interface, only: estimate_from_saved_orbitals
        use estimate_thomas_fermi_orbitals_interface, only: estimate_thomas_fermi_orbitals
        use estimate_hydrogenic_orbitals_interface, only: estimate_hydrogenic_orbitals
        use write_initial_wavefunction_summary_interface, only: write_initial_wavefunction_summary
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in) :: default_orbitals

        ! Local declarations
        integer  :: index(max_orbitals)
        integer :: j, k, nradial, nsubs, loc
        real(real64) :: con, fkk
        logical :: all, existed
        logical, save :: modify = .false.
        character(len=max(128, len(default_orbitals))) :: infile
        !-----------------------------------------------

        !
        !   Set the threshold for node counting
        !
        relative_extremum_cutoff = 0.05_real64
        !
        !   Set up the Thomas-Fermi potential
        !
        call build_thomas_fermi_potential
        !
        !   Allocate storage to the arrays that store the subshell
        !   radial wavefunction arrays; initialise these and all
        !   associated arrays
        !
        call alloc (radial_large_components,max_radial_points,orbital_count,"radial_large_components", &
            "generate_initial_radial_wavefunctions")
        call alloc (radial_small_components,max_radial_points,orbital_count,"radial_small_components", &
            "generate_initial_radial_wavefunctions")

        con = nuclear_charge/speed_of_light_au
        con = con*con

        do j = 1, orbital_count
            orbital_initialized(j) = .false.
            orbital_source_description(j) = " "
            radial_large_components(:radial_point_count,j) = 0.0_real64
            radial_small_components(:radial_point_count,j) = 0.0_real64

            k = abs(orbital_kappa(j))
            if (nuclear_model_parameter_count > 0) then
                orbital_origin_exponents(j) = dble(k)
            else if (nuclear_model_parameter_count == 0) then
                fkk = dble(k*k)
                if (fkk >= con) then
                    orbital_origin_exponents(j) = sqrt(fkk - con)
                else
                    write (stderr_unit, *) "generate_initial_radial_wavefunctions: Imaginary gamma parameter ", "for ", &
                        orbital_principal_n(j), orbital_symmetry_labels(j), " orbital;"
                    write (stderr_unit, *) "the point model for the nucleus ", &
                        "is inappropriate for Z > ", speed_of_light_au, "."
                    stop
                end if
            end if
        end do
        !
        !   Write out the complete list of subshell radial wave functions
        !
        call print_uninitialized_orbitals (all)
        !
        ! Direct to read radial functions till finish
        !
        do
            do while (.not.all)
                do
                    write (stderr_unit, *)
                    write (stderr_unit, *) "Read subshell radial wavefunctions. ", &
                        "Choose one below"
                    write (stderr_unit, *) "    1 -- GRASP92 File"
                    write (stderr_unit, *) "    2 -- Thomas-Fermi"
                    write (stderr_unit, *) "    3 -- Screened Hydrogenic"
                    write (stderr_unit, *) "    4 -- Screened Hydrogenic [custom Z]"

                    read (stdin_unit, *) nradial
                    if (nradial<1 .or. nradial>4) then
                        write (stderr_unit, *) nradial, "is not a valid choice, redo"
                    else
                        exit
                    end if
                end do

                if (nradial == 1) then
                    do
                        write (stderr_unit, *) 'Enter the file name (blank for ', default_orbitals, ')'
                        read (stdin_unit, "(A)") infile
                        if (len_trim(infile) == 0) infile = default_orbitals

                        inquire(file=infile, exist=existed)
                        if (.not.existed) then
                            write (stderr_unit, *) ' File "', infile(1:len_trim(infile)), &
                                '" does not exist, redo'
                        else
                            exit
                        end if
                    end do
                end if

                write (stderr_unit, *) "Enter the list of relativistic subshells:"
                open(newunit=calculation_log_unit, status="SCRATCH")
                call select_subshell_labels (index, nsubs)
                close(calculation_log_unit, status="DELETE")
                calculation_log_unit = -1

                if (nradial == 1) then
                    call estimate_from_saved_orbitals (index, nsubs, infile)
                else if (nradial == 2) then
                    call estimate_thomas_fermi_orbitals (index, nsubs)
                else if (nradial == 3) then
                    call estimate_hydrogenic_orbitals (index, nsubs, modify,nradial)
                else if (nradial == 4) then
                    call estimate_hydrogenic_orbitals (index, nsubs, modify,nradial)
                end if

                call print_uninitialized_orbitals (all)
            end do
            !
            ! All read. Let know, and allow modifying if non default
            !
            write (stderr_unit, *) "All required subshell radial wavefunctions ", &
                " have been estimated:"
            call write_initial_wavefunction_summary (stderr_unit)

            if (nondefault_options_flag == 0) then
                modify = .false.
            else
                write (stderr_unit, *) "Revise any of these estimates?"
                modify = read_yes_no()
            end if

            if (.not.modify) exit
            do
                write (stderr_unit, *) "Enter the list of subshells whose radial ", &
                    "wavefunctions are to be revised:"
                open(newunit=calculation_log_unit, status="SCRATCH")
                call select_subshell_labels (index, nsubs)
                close(calculation_log_unit, status="DELETE")
                calculation_log_unit = -1
                if (nsubs /= 0) exit
            end do
            do j = 1, nsubs
                loc = index(j)
                orbital_initialized(loc) = .false.
                radial_large_components(:radial_point_count,loc) = 0.0_real64
                radial_small_components(:radial_point_count,loc) = 0.0_real64
                orbital_source_description(loc) = " "
            end do
            all = .false.
        end do

        if (nondefault_options_flag /= 0) call write_initial_wavefunction_summary (summary_output_unit)

        return
    end subroutine generate_initial_radial_wavefunctions
end module generate_initial_radial_wavefunctions_interface
