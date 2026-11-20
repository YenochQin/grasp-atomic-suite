! Tabulate a point-Coulomb Dirac bound orbital, its binding energy, and its leading origin coefficient.
module dirac_coulomb_radial_orbital
    implicit none (type, external)

    ! Module interface
    private
    public :: dcbsrw
contains
    subroutine dcbsrw(principal_n, kappa, effective_nuclear_charge, binding_energy, large_component_origin_coefficient, &
        & large_component, small_component, last_radial_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters, only: max_radial_points
        use atomic_calculation_state ,          only: speed_of_light_au, accy
        use radial_grid_state, only: radial_grid, radial_point_count
        use radial_quadrature_workspace, only: first_series_coefficients => ta, second_series_coefficients => tb
        use complex_gamma, only: cgamma
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: principal_n
        integer, intent(in) :: kappa
        integer, intent(out) :: last_radial_index
        real(real64), intent(in) :: effective_nuclear_charge
        real(real64), intent(out) :: binding_energy
        real(real64), intent(out) :: large_component_origin_coefficient
        real(real64), intent(inout) :: large_component(max_radial_points)
        real(real64), intent(out) :: small_component(max_radial_points)

        ! Local declarations
        integer :: abs_kappa, series_index, radial_quantum_number, radial_quantum_factorial, sample_index, &
            & first_series_order, second_series_order
        real(real64) :: fine_structure_constant, principal_n_real, kappa_real, abs_kappa_real, radial_quantum_number_real, &
            & charge_times_alpha, origin_exponent, twice_exponent_plus_one, effective_principal_n, &
            relativistic_energy_ratio, gamma_argument_real, gamma_argument_imaginary, first_gamma_value, &
            & unused_gamma_imaginary, second_gamma_value, scale_factor, large_component_scale, small_component_scale, &
            & term_factorial, numerator_parameter, first_numerator_product, second_numerator_product&
            , denominator_parameter, denominator_product, inverse_series_denominator, effective_n_minus_kappa, scaled_radius, &
            & radius_power, first_radial_series, second_radial_series, radial_envelope, amplitude_cutoff
        if (principal_n <= 0) then
            write (*, '("DCBSRW:")')
            write (*, '(" Principal quantum number is ",1i3)') principal_n
            stop
        end if
        if (kappa == 0) then
            write (*, '("DCBSRW:")')
            write (*, '(" Angular quantum number is 0")')
            stop
        else if (kappa == principal_n) then
            write (*, '("DCBSRW:")')
            write (*, '(" Angular quantum number (",1i3,") is out of range for"," principal quantum number (",1i3,&
                &")")') kappa, principal_n
            stop
        else if (abs(kappa) > principal_n) then
            write (*, '("DCBSRW:")')
            write (*, '(" Angular quantum number (",1i3,") is out of range for"," principal quantum number (",1i3,&
                &")")') kappa, principal_n
            stop
        end if
        if (effective_nuclear_charge <= 0.0_real64) then
            write (*, '("DCBSRW:")')
            write (*, '(" Nuclear charge (",3p,1d16.7,") is too small")') effective_nuclear_charge
            stop
        else if (effective_nuclear_charge > speed_of_light_au) then
            write (*, '("DCBSRW:")')
            write (*, '(" Nuclear charge (",3p,1d16.7,") exceeds C (",1d16.7,")")') &
                effective_nuclear_charge, speed_of_light_au
            stop
        end if
        fine_structure_constant = 1.0_real64/speed_of_light_au
        principal_n_real = dble(principal_n)
        kappa_real = dble(kappa)
        abs_kappa = abs(kappa)
        abs_kappa_real = dble(abs_kappa)
        radial_quantum_number = principal_n - abs_kappa
        radial_quantum_number_real = dble(radial_quantum_number)
        charge_times_alpha = effective_nuclear_charge*fine_structure_constant
        origin_exponent = sqrt(abs_kappa_real*abs_kappa_real - charge_times_alpha*charge_times_alpha)
        twice_exponent_plus_one = origin_exponent + origin_exponent + 1.0_real64
        effective_principal_n = sqrt(principal_n_real*principal_n_real - &
            & 2.0_real64*radial_quantum_number_real*(abs_kappa_real - origin_exponent))
        relativistic_energy_ratio = 1.0_real64/sqrt(1.0_real64 + (charge_times_alpha/(origin_exponent + &
            & radial_quantum_number_real))**2)
        ! Return the positive binding energy relative to the continuum threshold.
        binding_energy = (1.0_real64 - relativistic_energy_ratio)*speed_of_light_au*speed_of_light_au
        radial_quantum_factorial = 1
        do sample_index = 1, radial_quantum_number
            radial_quantum_factorial = radial_quantum_factorial*sample_index
        end do
        gamma_argument_real = twice_exponent_plus_one + radial_quantum_number_real
        gamma_argument_imaginary = 0.0_real64
        call cgamma (gamma_argument_real, gamma_argument_imaginary, first_gamma_value, unused_gamma_imaginary)
        gamma_argument_real = twice_exponent_plus_one
        call cgamma (gamma_argument_real, gamma_argument_imaginary, second_gamma_value, unused_gamma_imaginary)
        scale_factor = -sqrt(first_gamma_value)/(second_gamma_value*sqrt(dble(radial_quantum_factorial)))* &
            & sqrt(effective_nuclear_charge/(2.0_real64*effective_principal_n*effective_principal_n*(&
            effective_principal_n - kappa_real)))
        if (kappa > 0) scale_factor = -scale_factor
        large_component_scale = scale_factor*sqrt(1.0_real64 + relativistic_energy_ratio)
        small_component_scale = scale_factor*sqrt(1.0_real64 - relativistic_energy_ratio)
        if (radial_quantum_number == 0) then
            first_series_order = 0
            second_series_order = 0
        else
            first_series_order = radial_quantum_number - 1
            second_series_order = radial_quantum_number
        end if
        scale_factor = 1.0_real64
        term_factorial = 1.0_real64
        numerator_parameter = -radial_quantum_number_real
        ! Generate the two terminating confluent hypergeometric series coefficients.
        first_numerator_product = numerator_parameter + 1.0_real64
        second_numerator_product = numerator_parameter
        denominator_parameter = twice_exponent_plus_one
        denominator_product = denominator_parameter
        series_index = 0
        do
            series_index = series_index + 1
            inverse_series_denominator = 1.0_real64/(term_factorial*denominator_product)
            if (series_index <= first_series_order) first_series_coefficients(series_index) = &
                & first_numerator_product*inverse_series_denominator
            if (series_index > second_series_order) exit
            second_series_coefficients(series_index) = second_numerator_product*inverse_series_denominator
            numerator_parameter = numerator_parameter + 1.0_real64
            first_numerator_product = first_numerator_product*(numerator_parameter + 1.0_real64)
            second_numerator_product = second_numerator_product*numerator_parameter
            denominator_parameter = denominator_parameter + 1.0_real64
            denominator_product = denominator_product*denominator_parameter
            scale_factor = scale_factor + 1.0_real64
            term_factorial = term_factorial*scale_factor
        end do
        large_component(1) = 0.0_real64
        small_component(1) = 0.0_real64
        scale_factor = (effective_nuclear_charge + effective_nuclear_charge)/effective_principal_n
        effective_n_minus_kappa = effective_principal_n - kappa_real
        do sample_index = 2, radial_point_count
            scaled_radius = scale_factor*radial_grid(sample_index)
            radius_power = scaled_radius
            series_index = 0
            first_radial_series = 1.0_real64
            second_radial_series = 1.0_real64
            do
                series_index = series_index + 1
                if (series_index <= first_series_order) first_radial_series = first_radial_series + &
                    & first_series_coefficients(series_index)*radius_power
                if (series_index > second_series_order) exit
                second_radial_series = second_radial_series + second_series_coefficients(series_index)*radius_power
                radius_power = radius_power*scaled_radius
            end do
            first_radial_series = radial_quantum_number_real*first_radial_series
            second_radial_series = effective_n_minus_kappa*second_radial_series
            radial_envelope = exp((-0.5_real64*scaled_radius))*scaled_radius**origin_exponent
            large_component(sample_index) = large_component_scale*radial_envelope*(first_radial_series - second_radial_series)
            small_component(sample_index) = small_component_scale*radial_envelope*(first_radial_series + second_radial_series)
        end do
        ! Trim from the outermost point using the large component and the hard radius cut.
        amplitude_cutoff = accy*0.1_real64
        last_radial_index = radial_point_count
        do while (abs(large_component(last_radial_index)) < amplitude_cutoff .or. &
            abs(radial_grid(last_radial_index)) > 1.0e50_real64)
            large_component(last_radial_index) = 0.0_real64
            small_component(last_radial_index) = 0.0_real64
            last_radial_index = last_radial_index - 1
        end do
        if (last_radial_index == radial_point_count) write (*, '(/,/,/," ***** Warning in SUBROUTINE DCBSRW *****"/,/&
            &," Radial grid of insufficient extent:"/," P(",1i4,") = ",1p,1d10.3,", Exceeds cutoff (",1d10.3,")")') &
            radial_point_count, large_component(radial_point_count), &
            & amplitude_cutoff
        large_component_origin_coefficient = large_component_scale*scale_factor**origin_exponent*(radial_quantum_number_real &
            & - effective_n_minus_kappa)
        return
        return
    end subroutine dcbsrw
end module dirac_coulomb_radial_orbital
