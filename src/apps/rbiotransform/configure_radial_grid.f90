module configure_radial_grid_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: configure_radial_grid
contains

    !***********************************************************************
    !                                                                      *
    subroutine configure_radial_grid
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: input_unit, output_unit
        !                                                                      *
        !   This subroutine sets the parameters controlling the radial grid    *
        !                                                                      *
        !   Last revision: June 1996                                           *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_radial_points
        use atomic_calculation_state,           only: speed_of_light_au, vacuum_speed_of_light_au, nuclear_charge, working_tolerance
        use calculation_input_defaults, only: nondefault_options_flag
        use radial_grid_state, only: radial_point_count, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count
        use yes_no_input, only: read_yes_no
        implicit none (type, external)

        ! Local declarations
        logical   :: yes
        !-----------------------------------------------
        !
        if (nuclear_model_parameter_count == 0) then
            radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
            transformed_grid_step = 0.5_real64**4
            radial_point_count = min(220,max_radial_points)
        else
            !cff     .. Should be Z-dependent
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

        return
    end subroutine configure_radial_grid
end module configure_radial_grid_interface
