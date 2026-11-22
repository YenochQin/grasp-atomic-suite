! Solve the inhomogeneous radial Dirac equations in the orbital tail region.
module integrate_radial_inward_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: integrate_radial_inward
contains

    ! Extend the method of C. Froese Fischer (Can. J. Phys. 41, 1895, 1963) to the inhomogeneous equations.
    ! Use the supplied large component at the matching point and find an outer point with a small tail.
    !
    ! The five-band finite-difference system M W = V couples alternating small and large components.
    ! Factor M by Crout elimination, solve L Z = V, then back-substitute in U W = Z.
    ! Last update: 10 Dec 1992
    subroutine integrate_radial_inward(orbital_index, matching_grid_index, large_component_values, small_component_values, &
        & last_radial_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use atomic_calculation_state,            only: working_tolerance
        use radial_grid_state, only: transformed_grid_step, radial_point_count, radial_derivative_over_radius
        use radial_dirac_workspace,            only: large_equation_coefficients, small_equation_coefficients, &
            large_radial_source, small_radial_source
        use orbital_csf_basis,            only: orbital_kappa
        use integer_format, only: format_integer
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index
        integer , intent(in) :: matching_grid_index
        integer , intent(out) :: last_radial_index
        real(real64) , intent(inout) :: large_component_values(max_radial_points)
        real(real64) , intent(inout) :: small_component_values(max_radial_points)

        ! Local declarations
        integer :: elimination_index, radial_index, text_length, first_zero_grid_index
        real(real64) :: small_component_pivots(max_radial_points), large_component_pivots(max_radial_points), &
            & large_equation_coupling(max_radial_points), small_component_rhs(max_radial_points), &
            large_component_rhs(max_radial_points)
        real(real64) :: tail_tolerance, half_step_times_kappa, radial_kappa_factor, one_plus_grid_kappa, &
            & one_minus_grid_kappa, next_one_plus_grid_kappa, next_one_minus_grid_kappa, &
            & large_from_current_small, current_tail_amplitude, small_from_previous_large, &
            large_from_previous_large, previous_tail_amplitude
        character :: number_text*5
        ! Set the tail tolerance and the scale for the kappa terms.
        tail_tolerance = 0.01_real64*working_tolerance
        half_step_times_kappa = 0.5_real64*transformed_grid_step*dble(orbital_kappa(orbital_index))
        ! Initialize counters
        elimination_index = 1
        radial_index = matching_grid_index

        radial_kappa_factor = half_step_times_kappa*radial_derivative_over_radius(radial_index)
        one_plus_grid_kappa = 1.0_real64 + radial_kappa_factor
        one_minus_grid_kappa = 1.0_real64 - radial_kappa_factor
        radial_kappa_factor = half_step_times_kappa*radial_derivative_over_radius(radial_index+1)
        next_one_plus_grid_kappa = 1.0_real64 + radial_kappa_factor
        next_one_minus_grid_kappa = 1.0_real64 - radial_kappa_factor
        ! Compute required elements of first two rows of  L  and  U
        small_component_pivots(elimination_index) = one_plus_grid_kappa
        large_from_current_small = -large_equation_coefficients(radial_index)/small_component_pivots(elimination_index)
        large_component_pivots(elimination_index) = (-next_one_plus_grid_kappa) + &
            & large_from_current_small*small_equation_coefficients(radial_index+1)
        large_equation_coupling(elimination_index) = next_one_minus_grid_kappa*large_from_current_small - &
            & large_equation_coefficients(radial_index+1)
        ! First elements of solution vector Z
        small_component_rhs(elimination_index) = (-small_radial_source(radial_index)) + small_equation_coefficients(radial_index)* &
            large_component_values(radial_index)
        large_component_rhs(elimination_index) = (-large_radial_source(radial_index)) - &
            & one_minus_grid_kappa*large_component_values(radial_index) - large_from_current_small* &
            small_component_rhs(elimination_index)
        current_tail_amplitude = abs(large_component_rhs(elimination_index)/large_component_pivots(elimination_index))
        do
            elimination_index = elimination_index + 1
            radial_index = radial_index + 1
            ! Failure if tables not long enough
            if (radial_index >= radial_point_count) then
                call format_integer (radial_point_count, number_text, text_length)
                write (output_unit, *) "IN: maximum tabulation point exceeds"
                write (output_unit, *) " dimensional limit (currently "//number_text(1:text_length)//");"
                write (output_unit, *) " radial wavefunction may indicate a"
                write (output_unit, *) " continuum state."
                error stop
            end if
            ! Compute required elements of remaining rows of  L  and  U
            one_plus_grid_kappa = next_one_plus_grid_kappa
            one_minus_grid_kappa = next_one_minus_grid_kappa
            radial_kappa_factor = half_step_times_kappa*radial_derivative_over_radius(radial_index+1)
            next_one_plus_grid_kappa = 1.0_real64 + radial_kappa_factor
            next_one_minus_grid_kappa = 1.0_real64 - radial_kappa_factor
            small_from_previous_large = -small_equation_coefficients(radial_index)/large_component_pivots(elimination_index-1)
            small_component_pivots(elimination_index) = one_plus_grid_kappa - &
                & small_from_previous_large*large_equation_coupling(elimination_index-1)
            large_from_previous_large = one_minus_grid_kappa/large_component_pivots(elimination_index-1)
            large_from_current_small = ((-large_equation_coefficients(radial_index)) - &
                & large_from_previous_large*large_equation_coupling(elimination_index-1)) &
                & /small_component_pivots(elimination_index)
            large_component_pivots(elimination_index) = (-next_one_plus_grid_kappa) + &
                & large_from_current_small*small_equation_coefficients(radial_index+1)
            large_equation_coupling(elimination_index) = (-large_equation_coefficients(radial_index+1)) + &
                & next_one_minus_grid_kappa*large_from_current_small
            ! Solution of  L*Z = V
            small_component_rhs(elimination_index) = (-small_radial_source(radial_index)) - &
                & small_from_previous_large*large_component_rhs(elimination_index-1)
            large_component_rhs(elimination_index) = (-large_radial_source(radial_index)) - &
                & large_from_previous_large*large_component_rhs(elimination_index-1) - &
                    & large_from_current_small*small_component_rhs(elimination_index)
            ! Test for outer boundary
            previous_tail_amplitude = current_tail_amplitude
            current_tail_amplitude = abs(large_component_rhs(elimination_index)/large_component_pivots(elimination_index))
            if (current_tail_amplitude + previous_tail_amplitude <= tail_tolerance .or. radial_index == max_radial_points-1 ) then
                last_radial_index = radial_index
                exit
            end if
        end do
        ! Reset counter
        elimination_index = elimination_index - 1
        ! Last two rows of solution of  U*W = Z ; evaluation of Q(J)
        small_component_values(radial_index) = 0.0_real64
        large_component_values(radial_index) = large_component_rhs(elimination_index)/large_component_pivots(elimination_index)
        small_component_values(radial_index-1) = &
            & (small_component_rhs(elimination_index) + &
                & small_equation_coefficients(radial_index)*large_component_values(radial_index)) &
            & /small_component_pivots(elimination_index)
        ! Solution of  U*W = Z
        do
            radial_index = radial_index - 1
            elimination_index = elimination_index - 1
            if (elimination_index <= 0) exit
            large_component_values(radial_index) = &
                & (large_component_rhs(elimination_index) &
                & - large_equation_coupling(elimination_index)*small_component_values(radial_index)) &
                & /large_component_pivots(elimination_index)
            small_component_values(radial_index-1) = &
                & (small_component_rhs(elimination_index) + &
                    & (1.0_real64-half_step_times_kappa*radial_derivative_over_radius(radial_index)) &
                & *small_component_values(radial_index) + &
                small_equation_coefficients(radial_index) &
                & *large_component_values(radial_index))/small_component_pivots(elimination_index)
        end do
        ! Complete tables with zeroes
        first_zero_grid_index = last_radial_index + 1
        large_component_values(first_zero_grid_index:radial_point_count) = 0.0_real64
        small_component_values(first_zero_grid_index:radial_point_count) = 0.0_real64
        return
    end subroutine integrate_radial_inward
end module integrate_radial_inward_interface
