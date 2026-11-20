! Initialize the first six radial samples for the homogeneous, inhomogeneous, or variational Dirac equation.
module radial_origin_solution
    implicit none (type, external)
    private
    public :: start
contains
      subroutine start(orbital_index, radial_equation_type, origin_slope_parameter, large_component_values, &
         & small_component_origin_coefficient, small_component_values)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use radial_difference_coefficients, only: cnc5c, six_point_weights => cnc6c
      use atomic_calculation_state,           only: speed_of_light_au, nuclear_charge, working_tolerance => accy
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: orbital_energy_parameters, orbital_kappa, orbital_origin_exponents
      use scf_radial_potentials,          only: orbital_y_potential => yp, large_exchange_source => xp, small_exchange_source => xq
      use scf_iteration_state,           only: direct_coefficient_count => ndcof, direct_orbital_indices => nda, &
         & direct_coefficients => da
      use radial_wavefunctions,          only: orbital_origin_parameters,radial_large_components,radial_small_components
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer, intent(in) :: radial_equation_type
      real(real64), intent(in) :: origin_slope_parameter
      real(real64), intent(out) :: small_component_origin_coefficient
      real(real64), intent(inout) :: large_component_values(max_radial_points)
      real(real64), intent(inout) :: small_component_values(max_radial_points)
      integer, parameter :: maximum_iterations = 36
      integer :: orbital_kappa_value, radial_or_coefficient_index, direct_orbital_index, iteration_count, start_grid_index
      real(real64) :: large_equation_coupling(6), small_equation_coupling(6), large_equation_source(6), small_equation_source(6)
      real(real64) :: large_component_estimates(2:6), small_component_estimates(2:6)
      real(real64) :: inverse_speed_of_light, charge_over_c, origin_exponent, kappa_real, one_minus_origin_exponent, &
         & one_minus_exponent_minus_kappa, one_minus_exponent_plus_kappa, large_origin_source_sum, &
         small_origin_source_sum, direct_orbital_origin_parameter, large_origin_value, small_origin_value, large_equation_sum, &
            & small_equation_sum, factor, c_squared, twice_c_squared, binding_energy, twice_c_squared_minus_energy&
         , radius, grid_derivative, radius_times_derivative, potential_times_derivative, unscaled_maximum_change, &
            & maximum_change, coefij, pi, qi, start_grid_radius
      inverse_speed_of_light = 1.0_real64/speed_of_light_au
      charge_over_c = nuclear_charge*inverse_speed_of_light
      origin_exponent = orbital_origin_exponents(orbital_index)
      orbital_kappa_value = orbital_kappa(orbital_index)
      kappa_real = dble(orbital_kappa_value)
      one_minus_origin_exponent = 1.0_real64 - origin_exponent
      one_minus_exponent_minus_kappa = one_minus_origin_exponent - kappa_real
      one_minus_exponent_plus_kappa = one_minus_origin_exponent + kappa_real
      ! The homogeneous and inhomogeneous modes use the point- or finite-nucleus origin coefficients.
      if (radial_equation_type==1 .or. radial_equation_type==2) then
         if (nuclear_model_parameter_count == 0) then
            large_component_values(1) = origin_slope_parameter
            if (orbital_kappa_value < 0) then
               small_component_values(1) = -origin_slope_parameter*charge_over_c/(origin_exponent - kappa_real)
            else
               small_component_values(1) = origin_slope_parameter*(origin_exponent + kappa_real)/charge_over_c
            end if
         else
            if (orbital_kappa_value < 0) then
               large_component_values(1) = origin_slope_parameter
               small_component_values(1) = 0.0_real64
            else
               large_component_values(1) = 0.0_real64
               small_component_values(1) = origin_slope_parameter*(origin_exponent + kappa_real)/charge_over_c
            end if
         end if
         if (radial_equation_type == 1) then
            large_equation_source = 0.0_real64
            small_equation_source = 0.0_real64
         else
            large_origin_source_sum = 0.0_real64
            small_origin_source_sum = 0.0_real64
            do radial_or_coefficient_index = 1, direct_coefficient_count
               direct_orbital_index = direct_orbital_indices(radial_or_coefficient_index)
               direct_orbital_origin_parameter = orbital_origin_parameters(direct_orbital_index)
               if (nuclear_model_parameter_count == 0) then
                  large_origin_value = direct_orbital_origin_parameter
                  if (orbital_kappa_value < 0) then
                     small_origin_value = -direct_orbital_origin_parameter*charge_over_c/(origin_exponent - kappa_real)
                  else
                     small_origin_value = direct_orbital_origin_parameter*(origin_exponent + kappa_real)/charge_over_c
                  end if
                  large_equation_sum = charge_over_c*small_origin_value
                  small_equation_sum = -charge_over_c*large_origin_value
               else
                  if (orbital_kappa_value < 0) then
                     large_origin_value = direct_orbital_origin_parameter
                     small_origin_value = 0.0_real64
                  else
                     large_origin_value = 0.0_real64
                     small_origin_value = direct_orbital_origin_parameter*(origin_exponent + kappa_real)/charge_over_c
                  end if
                  large_equation_sum = 0.0_real64
                  small_equation_sum = 0.0_real64
               end if
               factor = direct_coefficients(radial_or_coefficient_index)
               large_origin_source_sum = large_origin_source_sum + factor*(large_equation_sum + &
                  & one_minus_exponent_minus_kappa*large_origin_value)
               small_origin_source_sum = small_origin_source_sum + factor*(small_equation_sum + &
                  & one_minus_exponent_plus_kappa*small_origin_value)
            end do
            factor = radial_grid_derivative(1)
            large_equation_source(1) = factor*large_origin_source_sum
            small_equation_source(1) = factor*small_origin_source_sum
            do radial_or_coefficient_index = 2, 6
               factor = -radial_grid_derivative(radial_or_coefficient_index)*radial_grid(radial_or_coefficient_index)**(- &
                  & origin_exponent)
               large_equation_source(radial_or_coefficient_index) = factor*large_exchange_source(radial_or_coefficient_index)
               small_equation_source(radial_or_coefficient_index) = factor*small_exchange_source(radial_or_coefficient_index)
            end do
         end if
      ! The variational mode starts at zero and builds its source from the current orbital.
      else if (radial_equation_type == 3) then
         large_component_values(1) = 0.0_real64
         small_component_values(1) = 0.0_real64
         large_equation_source(1) = 0.0_real64
         small_equation_source(1) = 0.0_real64
         do radial_or_coefficient_index = 2, 6
            factor = inverse_speed_of_light*radial_grid_derivative(radial_or_coefficient_index)* &
               & radial_grid(radial_or_coefficient_index)**one_minus_origin_exponent
            large_equation_source(radial_or_coefficient_index) = &
               & -factor*radial_small_components(radial_or_coefficient_index,orbital_index)
            small_equation_source(radial_or_coefficient_index) = factor*radial_large_components(radial_or_coefficient_index, &
               & orbital_index)
         end do
      end if
      small_component_origin_coefficient = small_component_values(1)
      c_squared = speed_of_light_au*speed_of_light_au
      twice_c_squared = c_squared + c_squared
      binding_energy = orbital_energy_parameters(orbital_index)
      twice_c_squared_minus_energy = twice_c_squared - binding_energy
      do radial_or_coefficient_index = 1, 6
         radius = radial_grid(radial_or_coefficient_index)
         grid_derivative = radial_grid_derivative(radial_or_coefficient_index)
         radius_times_derivative = radius*grid_derivative
         potential_times_derivative = orbital_y_potential(radial_or_coefficient_index)*grid_derivative
         large_equation_coupling(radial_or_coefficient_index) = &
            & -inverse_speed_of_light*(twice_c_squared_minus_energy*radius_times_derivative + potential_times_derivative)
         small_equation_coupling(radial_or_coefficient_index) = &
            & -inverse_speed_of_light*(binding_energy*radius_times_derivative - potential_times_derivative)
      end do
      iteration_count = 0
      large_origin_value = large_component_values(1)
      small_origin_value = small_component_values(1)
      unscaled_maximum_change = max(abs(large_origin_value),abs(small_origin_value))
      large_component_values(2:6) = large_origin_value
      small_component_values(2:6) = small_origin_value
      ! Measure convergence in the scaled-back radial components at the edge of the starting interval.
      factor = radial_grid(6)**origin_exponent
    7 continue
      iteration_count = iteration_count + 1
      maximum_change = 0.0_real64
      do start_grid_index = 2, 6
         large_equation_sum = sum(six_point_weights(:, &
            & start_grid_index)*(one_minus_exponent_minus_kappa*radial_grid_derivative(:6)*large_component_values(:6)- &
            & large_equation_coupling*small_component_values(:6)+large_equation_source))
         small_equation_sum = sum(six_point_weights(:, &
            & start_grid_index)*(one_minus_exponent_plus_kappa*radial_grid_derivative(:6)*small_component_values(:6)- &
            & small_equation_coupling*large_component_values(:6)+small_equation_source))
         start_grid_radius = radial_grid(start_grid_index)
         large_equation_sum = large_equation_sum/start_grid_radius
         small_equation_sum = small_equation_sum/start_grid_radius
         large_component_estimates(start_grid_index) = large_equation_sum
         small_component_estimates(start_grid_index) = small_equation_sum
         maximum_change = max(maximum_change,abs(large_equation_sum - large_component_values(start_grid_index)))
         maximum_change = max(maximum_change,abs(small_equation_sum - small_component_values(start_grid_index)))
      end do
      large_component_values(2:6) = large_component_estimates(:6)
      small_component_values(2:6) = small_component_estimates(:6)
      unscaled_maximum_change = maximum_change
      maximum_change = maximum_change*factor
      if (maximum_change > working_tolerance) then
         if (iteration_count < maximum_iterations) then
            goto 7
         end if
      end if
      if (maximum_change > working_tolerance) then
         large_component_values(2:6) = large_origin_value
         small_component_values(2:6) = small_origin_value
      end if
      large_component_values(1) = 0.0_real64
      small_component_values(1) = 0.0_real64
      do radial_or_coefficient_index = 2, 6
         factor = radial_grid(radial_or_coefficient_index)**origin_exponent
         large_component_values(radial_or_coefficient_index) = factor*large_component_values(radial_or_coefficient_index)
         small_component_values(radial_or_coefficient_index) = factor*small_component_values(radial_or_coefficient_index)
      end do
      return
      return
      end subroutine start
end module radial_origin_solution
