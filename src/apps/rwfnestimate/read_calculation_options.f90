! Read the settings that select CSF interactions for angular-coefficient generation.
module read_calculation_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_calculation_options
contains

    ! Interactively determines data useful for generating estimates of
    ! the subshell radial wavefunctions.
    !
    ! Written by Farid A. Parpia. Last revision: 21 Dec 1992
    subroutine read_calculation_options
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use standard_io_units, only: stderr_unit
        use atomic_calculation_state, only: working_tolerance, max_radial_points, nuclear_charge, speed_of_light_au, &
            vacuum_speed_of_light_au
        use calculation_input_defaults, only: nondefault_options_flag
        use radial_grid_state, only: radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use yes_no_input, only: read_yes_no
        use isotope_data_input, only: read_isotope_data
        use radial_stencil_weights, only: initialize_radial_stencil_weights
        use radial_grid_setup, only: initialize_radial_grid
        use nuclear_coulomb_potential, only: calculate_nuclear_coulomb_potential
        implicit none (type, external)

        ! Local declarations
        !
        !   G l o b a l   P a r a m e t e r s
        logical :: yes
        ! Read the nuclear parameters from isodata.
        call read_isotope_data ("isodata")
        ! Set defaults
        speed_of_light_au = vacuum_speed_of_light_au
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
            write (stderr_unit, *) "Change the default speed of light ", &
                "or radial grid parameters?"
            yes = read_yes_no()
            if (yes) then
                ! Modify the speed of light
                write (stderr_unit, *) "The physical speed of light in ", &
                    "atomic units is", vacuum_speed_of_light_au, ";"
                write (stderr_unit, *) "revise this value?"
                yes = read_yes_no()
                if (yes) then
                    write (stderr_unit, *) "Enter the revised value:"
                    read (input_unit, *) speed_of_light_au
                end if
                ! Modify the parameters controlling the radial grid
                write (stderr_unit, *) "The default radial grid parameters ", &
                    "for this case are:"
                write (stderr_unit, *) " RNT = ", radial_grid_scale, ";"
                write (stderr_unit, *) " H = ", transformed_grid_step, ";"
                write (stderr_unit, *) " HP = ", asymptotic_grid_step, ";"
                write (stderr_unit, *) " N = ", radial_point_count, ";"
                write (stderr_unit, *) "revise these values?"
                yes = read_yes_no()
                if (yes) then
                    write (stderr_unit, *) "Enter RNT:"
                    read (input_unit, *) radial_grid_scale
                    write (stderr_unit, *) "Enter H:"
                    read (input_unit, *) transformed_grid_step
                    write (stderr_unit, *) "Enter HP:"
                    read (input_unit, *) asymptotic_grid_step
                    write (stderr_unit, *) "Enter N:"
                    read (input_unit, *) radial_point_count
                end if
            end if
        end if
        ! ACCY is an estimate of the accuracy of the numerical procedures
        working_tolerance = transformed_grid_step**6
        ! Initialize the radial integration and derivative stencils.
        call initialize_radial_stencil_weights
        ! Build the radial grid and its derivative arrays.
        call initialize_radial_grid
        ! Generate $- r \times V_ (r)$
        call calculate_nuclear_coulomb_potential
        return
    end subroutine read_calculation_options
end module read_calculation_options_interface
