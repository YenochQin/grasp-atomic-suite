! Construct the transformed radial grid and its Jacobian, including the workspace padding.
module radial_grid_setup
    implicit none (type, external)

    ! Module interface
    private
    public :: initialize_radial_grid
contains
    subroutine initialize_radial_grid
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use debug_output_settings, only: rmcdhf_debug_flags
        use atomic_calculation_state, only: machine_precision
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_derivative, radial_derivative_over_radius, &
            radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        implicit none (type, external)

        ! Local declarations
        integer :: grid_points_with_padding, radial_index, half_grid_point_count, debug_row_count, debug_row_index, &
            & left_grid_index, right_grid_index
        real(real64) :: exp_grid_step, exp_grid_coordinate, exp_grid_coordinate_minus_one, newton_relative_tolerance, &
            & grid_step_ratio, previous_radius_estimate, radius_estimate, transformed_coordinate, shifted_radius_estimate, &
            & grid_equation_residual, &
            inverse_grid_equation_derivative, radius_correction
        ! Initialize the extra points used by the radial stencils as well as the active grid.
        grid_points_with_padding = radial_point_count + 10
        radial_grid(1) = 0.0_real64
        radial_derivative_over_radius(1) = 0.0_real64
        ! Use the exponential transform when the asymptotic step is zero. Otherwise solve the mixed grid equation.
        if (asymptotic_grid_step == 0.0_real64) then
            radial_grid_derivative(1) = radial_grid_scale
            exp_grid_step = exp(transformed_grid_step)
            exp_grid_coordinate = 1.0_real64
            do radial_index = 2, grid_points_with_padding
                exp_grid_coordinate = exp_grid_step*exp_grid_coordinate
                exp_grid_coordinate_minus_one = exp_grid_coordinate - 1.0_real64
                radial_grid(radial_index) = radial_grid_scale*exp_grid_coordinate_minus_one
                radial_grid_derivative(radial_index) = radial_grid_scale*exp_grid_coordinate
                radial_derivative_over_radius(radial_index) = exp_grid_coordinate/exp_grid_coordinate_minus_one
            end do
        else
            ! Solve t = log(1+r/r0) + (h/hp)*r. The inverse slope is the stored dr/dt.
            newton_relative_tolerance = 1.0e3_real64*machine_precision
            grid_step_ratio = transformed_grid_step/asymptotic_grid_step
            radial_grid_derivative(1) = radial_grid_scale/(grid_step_ratio*radial_grid_scale + 1.0_real64)
            previous_radius_estimate = 0.0_real64
            radius_estimate = 0.0_real64
            do radial_index = 2, grid_points_with_padding
                transformed_coordinate = transformed_grid_step*dble(radial_index - 1)
                do
                    shifted_radius_estimate = radius_estimate + radial_grid_scale
                    grid_equation_residual = log(shifted_radius_estimate/radial_grid_scale) + grid_step_ratio*radius_estimate &
                        & - transformed_coordinate
                    inverse_grid_equation_derivative = shifted_radius_estimate / &
                        (grid_step_ratio*shifted_radius_estimate + 1.0_real64)
                    radius_correction = -grid_equation_residual*inverse_grid_equation_derivative
                    radius_estimate = previous_radius_estimate + radius_correction
                    if (abs(radius_correction/radius_estimate) < newton_relative_tolerance) exit
                    previous_radius_estimate = radius_estimate
                end do
                radial_grid(radial_index) = radius_estimate
                shifted_radius_estimate = radius_estimate + radial_grid_scale
                inverse_grid_equation_derivative = shifted_radius_estimate/(grid_step_ratio*shifted_radius_estimate + 1.0_real64)
                radial_grid_derivative(radial_index) = inverse_grid_equation_derivative
                radial_derivative_over_radius(radial_index) = inverse_grid_equation_derivative/radius_estimate
            end do
        end if
        if (rmcdhf_debug_flags(1)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,"From SUBROUTINE RADGRD:"/,2(" -------- r -------- -------- r'' -------",&
                &" ------- r''/r ------"))')
            half_grid_point_count = radial_point_count/2
            if (2*half_grid_point_count == radial_point_count) then
                debug_row_count = half_grid_point_count
            else
                debug_row_count = half_grid_point_count + 1
            end if
            do debug_row_index = 1, debug_row_count
                left_grid_index = debug_row_index
                right_grid_index = left_grid_index + debug_row_count
                if (right_grid_index <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d19.12))') &
                        radial_grid(left_grid_index), radial_grid_derivative(left_grid_index), &
                        & radial_derivative_over_radius(left_grid_index), radial_grid(right_grid_index), &
                        & radial_grid_derivative(right_grid_index), &
                        radial_derivative_over_radius(right_grid_index)
                else if (left_grid_index <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d19.12))') &
                        radial_grid(left_grid_index), radial_grid_derivative(left_grid_index), &
                        & radial_derivative_over_radius(left_grid_index)
                end if
            end do
        end if
        return
        return
    end subroutine initialize_radial_grid
end module radial_grid_setup
