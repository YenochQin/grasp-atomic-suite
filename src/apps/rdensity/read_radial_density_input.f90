module read_radial_density_input_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_radial_density_input
contains

    ! Interactively determines the data governing the SMS problem.
    !
    ! Written by Farid A. Parpia. Last revision: 15 Dec 1992
    subroutine read_radial_density_input(name)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use calculation_state_files, only: orbital_input_path
        use grasp_parameters,    only: max_radial_points
        use hamiltonian_correction_flags, only: use_first_order_csf_perturbation
        use atomic_calculation_state, only: working_tolerance, max_radial_points, nuclear_charge, speed_of_light_au, &
            vacuum_speed_of_light_au
        use calculation_input_defaults, only: nondefault_options_flag
        use first_order_csf_cutoffs, only: first_order_csf_cutoff
        use radial_grid_state, only: radial_point_count, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use nuclear_potential_state, only: max_radial_points
        use yes_no_input, only: read_yes_no
        use isotope_data_input, only: read_isotope_data
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_grid_setup, only: initialize_radial_grid
        use nuclear_coulomb_potential, only: calculate_nuclear_coulomb_potential
        use radial_wavefunction_input, only: read_radial_wavefunction_file
        implicit none (type, external)

        ! Arguments
        character, intent(in) :: name*24

        ! Local declarations
        logical :: yes
        ! Read the nuclear parameters from isodata.
        call read_isotope_data ("isodata")
        ! Determine the physical effects specifications

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
        if (nondefault_options_flag /= 0) then
            write (output_unit, *) "Treat contributions of some CSFs"
            write (output_unit, *) " as first-order perturbations?"
            yes = read_yes_no()
            if (yes) then
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
        ! Choose the radial grid parameters.
        if (nuclear_model_parameter_count == 0) then
            radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
            transformed_grid_step = 0.5_real64**4
            radial_point_count = min(220,max_radial_points)
        else
            ! The screening estimate could be refined to depend on nuclear charge Z.
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
        ! ACCY is an estimate of the accuracy of the numerical procedures
        working_tolerance = transformed_grid_step**6
        ! Initialize the radial integration and derivative stencils.
        call initialize_radial_stencil_weights
        ! Build the radial grid and its derivative arrays.
        call initialize_radial_grid
        ! Generate $- r \times V_nuc (r)$
        !
        ! Load the radial wavefunctions
        call read_radial_wavefunction_file (orbital_input_path(name))
        return
    end subroutine read_radial_density_input
end module read_radial_density_input_interface
