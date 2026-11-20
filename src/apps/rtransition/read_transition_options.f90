module read_transition_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_transition_options
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_transition_options(isotope_file_path,name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
        !                                                                      *
        !   Interactively determines the data governing the transition prob-   *
        !   lem.                                                               *
        !                                                                      *
        !   Written by Per Jonsson                Last revision: June 1996     *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: max_radial_points
        use atomic_calculation_state, only: working_tolerance, max_radial_points, nuclear_charge, speed_of_light_au, &
            vacuum_speed_of_light_au
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation
        use calculation_input_defaults, only: nondefault_options_flag
        use first_order_csf_cutoffs, only: first_order_csf_cutoff
        use radial_grid_state, only: transformed_grid_step, asymptotic_grid_step, radial_point_count, radial_grid_scale
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use transition_matrix_workspace, only: transition_option_flags
        use yes_no_input, only: read_yes_no
        use isotope_data_input, only: read_isotope_data

        ! Procedure interfaces
        use read_multipole_options_interface, only: read_multipole_options

        ! Module dependencies
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_grid_setup, only: initialize_radial_grid

        ! Procedure interfaces
        use read_initial_radial_wavefunctions_interface, only: read_initial_radial_wavefunctions
        use read_final_radial_wavefunctions_interface, only: read_final_radial_wavefunctions
        use calculate_biorthogonal_overlap_interface, only: calculate_biorthogonal_overlap
        implicit none (type, external)

        ! Arguments
        character  :: name(2)*24
        character(len=*) :: isotope_file_path

        ! Local declarations
        logical :: yes
        character :: cunits*4
        !-----------------------------------------------
        !
        !   Open, check, load data from, and close the  .iso  file
        !
        !GG      CALL SETISO ('isodata')
        call read_isotope_data (isotope_file_path)
        !
        !   Determine the physical effects specifications
        !
        if (nondefault_options_flag /= 0) then
            write (output_unit, *) "The physical speed of light in"
            write (output_unit, *) " atomic units is", vacuum_speed_of_light_au, ";"
            write (output_unit, *) " revise this value?"
            yes = read_yes_no()
            if (yes) then
                write (output_unit, *) "Enter the revised value:"
                read (input_unit, *) speed_of_light_au
            else
                speed_of_light_au = vacuum_speed_of_light_au
            end if
        else
            speed_of_light_au = vacuum_speed_of_light_au
        end if
        !
        use_first_order_csf_perturbation = .false.
        first_order_csf_cutoff = 0
        !
        !   Determine the multipolarity and parity of the transitions
        !
        call read_multipole_options
        !
        !   Determine the units for the printout and other options
        !
        transition_option_flags = .false.
        !
        if (nondefault_options_flag /= 0) then
            write (output_unit, *) "Which units are to be used to"
            write (output_unit, *) " express the transition energies?"
            write (output_unit, *) "    A    : Angstrom:"
            write (output_unit, *) "    eV   : electron volts;"
            write (output_unit, *) "    Hart : Hartree atomic units;"
            write (output_unit, *) "    Hz   : Hertz;"
            write (output_unit, *) "    Kays : Kaysers [cm**(-1)];"
            do
                read (*, "(A)") cunits
                if (cunits(1:1) == "A") then
                    transition_option_flags(1) = .true.
                else if (cunits(1:2) == "eV") then
                    transition_option_flags(2) = .true.
                else if (cunits(1:4) == "Hart") then
                    transition_option_flags(3) = .true.
                else if (cunits(1:2) == "Hz") then
                    transition_option_flags(4) = .true.
                else if (cunits(1:4) == "Kays") then
                    transition_option_flags(5) = .true.
                else
                    write (output_unit, *) "GETOSD: Unable to interpret string;"
                    write (output_unit, *) " reenter ..."
                    cycle
                end if
                exit
            end do
        else
            transition_option_flags(5) = .true.
        end if
        !
        !      WRITE (6, *) 'Sort transitions by energy?'
        !      YES = GETYN()
        !      IF (YES) LTC(6) = .TRUE.
        !
        if (nondefault_options_flag /= 0) then
            write (output_unit, *) "Einstein A and B coefficients are"
            write (output_unit, *) " printed in SI units; use Hartree"
            write (output_unit, *) " atomic units instead?"
            yes = read_yes_no()
            if (yes) transition_option_flags(7) = .true.
        else
            transition_option_flags(7) = .false.
        end if
        !
        !   Determine the parameters controlling the radial grid
        !
        if (nuclear_model_parameter_count == 0) then
            radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
            transformed_grid_step = 0.5_real64**4
            radial_point_count = min(220,max_radial_points)
        else
            !CFF     .. should be Z-dependent
            radial_grid_scale = 2.0e-6_real64/nuclear_charge
            transformed_grid_step = 5.0e-2_real64
            radial_point_count = max_radial_points
        end if
        asymptotic_grid_step = 0.0_real64
        if (nondefault_options_flag /= 0) then
            write (output_unit, *) "The default radial grid parameters"
            write (output_unit, *) " for this case are:"
            write (output_unit, *) " RNT = ", radial_grid_scale, ";"
            write (output_unit, *) " H = ", transformed_grid_step, ";"
            write (output_unit, *) " HP = ", asymptotic_grid_step, ";"
            write (output_unit, *) " N = ", radial_point_count, ";"
            write (output_unit, *) " revise these values?"
            yes = read_yes_no()
            if (yes) then
                write (output_unit, *) "Enter RNT:"
                read (input_unit, *) radial_grid_scale
                write (output_unit, *) "Enter H:"
                read (input_unit, *) transformed_grid_step
                write (output_unit, *) "Enter HP:"
                read (input_unit, *) asymptotic_grid_step
                write (output_unit, *) "Enter N:"
                read (input_unit, *) radial_point_count
            end if
        end if
        !
        !   ACCY is an estimate of the accuracy of the numerical procedures
        !
        working_tolerance = transformed_grid_step**6
        !
        !   Set up the coefficients for the numerical proceduRes
        !
        call initialize_radial_stencil_weights
        !
        !   Generate the radial grid and all associated arrays
        !
        call initialize_radial_grid
        !
        !   Load the initial state radial wavefunctions.
        !
        call read_initial_radial_wavefunctions (name(1))
        !
        !   Load the final state radial wavefunctions.
        !
        call read_final_radial_wavefunctions (name(2))

        !   Construct the radial overlap matrix.
        !
        call calculate_biorthogonal_overlap

        return
    end subroutine read_transition_options
end module read_transition_options_interface
