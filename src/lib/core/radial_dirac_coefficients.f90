! Prepare Dirac integration coefficients on the radial grid and choose an inward/outward matching point.
module radial_dirac_coefficients
    implicit none (type, external)
    private
    public :: setpot
contains
      subroutine setpot(orbital_index, matching_grid_index)
      use, intrinsic :: iso_fortran_env, only: real64
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use atomic_calculation_state, only: working_tolerance => accy, max_radial_points, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_dirac_workspace,           only: large_equation_coefficients => tf, small_equation_coefficients => tg
      use orbital_csf_basis, only: orbital_energy_parameters
      use scf_radial_potentials,          only: orbital_y_potential => yp
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer, intent(out) :: matching_grid_index
      integer :: radial_index
      real(real64) :: negative_step_over_twice_c, binding_energy, twice_c_squared_minus_energy, grid_derivative, &
         & scaled_orbital_potential
      logical :: searching_for_match
      searching_for_match = .true.
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
         if (.not.searching_for_match) cycle
         if (abs(sign(1.0_real64,small_equation_coefficients(radial_index))+sign(1.0_real64, &
            & small_equation_coefficients(1))) >= working_tolerance) cycle
         matching_grid_index = radial_index
         searching_for_match = .false.
      end do
      if (matching_grid_index == 0) then
         write (stderr_unit, *) "SETPOT: Join set to NNNP/2 = ", max_radial_points/2
         matching_grid_index = max_radial_points/2
      end if
      return
      end subroutine setpot
end module radial_dirac_coefficients
