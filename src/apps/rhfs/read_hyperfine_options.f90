module read_hyperfine_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_hyperfine_options
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_hyperfine_options(name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
        !                                                                      *
        !   Interactively determines the data governing the HFS problem.       *
        !                                                                      *
        !   Call(s) to: [LIB92]: NUCPOT, RADGRD, SETQIC.                       *
        !               [RCI92]: SETISO, SETRWF.                               *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use calculation_state_files, only: orbital_input_path
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation
        use atomic_calculation_state, only: working_tolerance, max_radial_points, nuclear_charge, speed_of_light_au, &
            vacuum_speed_of_light_au
        use calculation_input_defaults, only: nondefault_options_flag
        use first_order_csf_cutoffs, only: first_order_csf_cutoff
        use radial_grid_state,          only: radial_grid_scale, transformed_grid_step, radial_point_count, asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use yes_no_input, only: read_yes_no
        use isotope_data_input, only: read_isotope_data
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_grid_setup, only: initialize_radial_grid
        use nuclear_coulomb_potential, only: calculate_nuclear_coulomb_potential
        use radial_wavefunction_input, only: read_radial_wavefunction_file
        implicit none (type, external)

        ! Arguments
        character(len=24), intent(in) :: name

        ! Local declarations
        logical            :: prompt_answer
        character(len=128) :: isotope_file_path
        !-----------------------------------------------
        !
        !   Open, check, load data from, and close the  .iso  file
        !
        isotope_file_path = "isodata"
        call read_isotope_data (isotope_file_path)
        !
        !   Determine the physical effects specifications
        !
        if (nondefault_options_flag /= 0) then
            write (output_unit, *) "The physical speed of light in"
            write (output_unit, *) " atomic units is", vacuum_speed_of_light_au, ";"
            write (output_unit, *) " revise this value?"
            prompt_answer = read_yes_no()
            if (prompt_answer) then
                write (output_unit, *) "Enter the revised value:"
                read (input_unit, *) speed_of_light_au
            else
                speed_of_light_au = vacuum_speed_of_light_au
            end if
        else
            speed_of_light_au = vacuum_speed_of_light_au
        end if
        !
        if (nondefault_options_flag /= 0) then
            !
            write (output_unit, *) "Treat contributions of some CSFs"
            write (output_unit, *) " as first-order perturbations?"
            prompt_answer = read_yes_no()
            if (prompt_answer) then
                use_first_order_csf_perturbation = .true.
                write (output_unit, *) "The contribution of CSFs"
                write (output_unit, *) " 1 -- ICCUT will be treated"
                write (output_unit, *) " variationally; the remainder"
                write (output_unit, *) " perturbatively; enter ICCUT:"
                read (input_unit, *) first_order_csf_cutoff
            else
                use_first_order_csf_perturbation = .false.
                first_order_csf_cutoff = 0
            end if
        else
            use_first_order_csf_perturbation = .false.
            first_order_csf_cutoff = 0
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
            prompt_answer = read_yes_no()
            if (prompt_answer) then
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
        !   Set up the coefficients for the numerical procedures
        !
        call initialize_radial_stencil_weights
        !
        !   Generate the radial grid and all associated arrays
        !
        call initialize_radial_grid
        !
        !   Generate $- r \times V_nuc (r)$
        !
        call calculate_nuclear_coulomb_potential
        !
        !   Load the radial wavefunctions
        !
        !      CALL SETRWFA(NAME)

        call read_radial_wavefunction_file (orbital_input_path(name))

        !
        return
    end subroutine read_hyperfine_options
end module read_hyperfine_options_interface
