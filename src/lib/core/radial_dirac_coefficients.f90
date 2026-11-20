! Prepare Dirac integration coefficients on the radial grid and choose an inward/outward matching point.
module radial_dirac_coefficients
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_radial_dirac_coefficients
contains
    subroutine calculate_radial_dirac_coefficients(orbital_index, matching_grid_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use standard_io_units, only: stderr_unit
        use atomic_calculation_state, only: working_tolerance, max_radial_points, speed_of_light_au
        use radial_grid_state, only: radial_point_count, radial_grid_derivative, radial_derivative_over_radius, &
            transformed_grid_step
        use radial_dirac_workspace,           only: large_equation_coefficients, small_equation_coefficients
        use orbital_csf_basis, only: orbital_energy_parameters
        use scf_radial_potentials,          only: orbital_y_potential
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        integer, intent(out) :: matching_grid_index

        ! Local declarations
        integer :: radial_index
        real(real64) :: negative_step_over_twice_c, binding_energy, twice_c_squared_minus_energy, grid_derivative, &
            & scaled_orbital_potential
        logical :: searching_for_match
        searching_for_match = .true.
        ! Store step-scaled off-diagonal coefficients for integration in the transformed coordinate.
        negative_step_over_twice_c = -transformed_grid_step/(2.0_real64*speed_of_light_au)
        binding_energy = orbital_energy_parameters(orbital_index)
        twice_c_squared_minus_energy = 2.0_real64*speed_of_light_au*speed_of_light_au - binding_energy
        matching_grid_index = 0
        large_equation_coefficients(1) = 0.0_real64
        small_equation_coefficients(1) = 0.0_real64
        do radial_index = 2, radial_point_count
            grid_derivative = radial_grid_derivative(radial_index)
            scaled_orbital_potential = orbital_y_potential(radial_index)*radial_derivative_over_radius(radial_index)
            large_equation_coefficients(radial_index) = &
                & negative_step_over_twice_c*(twice_c_squared_minus_energy*grid_derivative + scaled_orbital_potential)
            small_equation_coefficients(radial_index) = negative_step_over_twice_c*(binding_energy*grid_derivative - &
                & scaled_orbital_potential)
            ! Join the inward and outward solutions at the first sign change relative to the origin.
            if (.not.searching_for_match) cycle
            if (abs(sign(1.0_real64,small_equation_coefficients(radial_index))+sign(1.0_real64, &
                & small_equation_coefficients(1))) >= working_tolerance) cycle
            matching_grid_index = radial_index
            searching_for_match = .false.
        end do
        if (matching_grid_index == 0) then
            write (stderr_unit, *) "SETPOT: Join set to max_radial_points/2 = ", max_radial_points/2
            matching_grid_index = max_radial_points/2
        end if
        return
    end subroutine calculate_radial_dirac_coefficients
end module radial_dirac_coefficients
