! Tabulate -r*V(r) for a point nucleus or a two-parameter Fermi charge distribution.
module nuclear_coulomb_potential
    implicit none (type, external)
    private
    public :: nucpot
contains
      subroutine nucpot
      use, intrinsic :: iso_fortran_env, only: real64
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state,           only: nuclear_charge, pi
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_potential_state,          only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
      use alternating_exponential_series, only: es
      implicit none (type, external)
      integer :: radial_index, points_per_column, table_row_count, table_row, first_column_index, second_column_index, &
         & third_column_index
      real(real64) :: half_density_radius, diffuseness, diffuseness_over_radius, twice_a_over_c, a_over_c_squared, &
         & three_a_over_c_squared, a_over_c_cubed, radius_over_diffuseness, pi_squared, half_pi_squared_a_over_c_squared, &
         interior_linear_factor, origin_inverse_square_sum, origin_inverse_cube_sum, six_a_over_c_cubed, &
            & origin_series_correction, fermi_normalization, normalized_nuclear_charge, radius, radius_from_half_density, &
            & scaled_radius_offset, radius_over_half_density, &
         radius_inverse_square_sum, radius_inverse_cube_sum
      logical :: found_first_exterior_point
      if (nuclear_model_parameter_count == 0) then
         negative_radius_times_potential(:radial_point_count) = nuclear_charge
      else if (nuclear_model_parameter_count == 2) then
         half_density_radius = nuclear_model_parameters(1)
         diffuseness = nuclear_model_parameters(2)
         diffuseness_over_radius = diffuseness/half_density_radius
         twice_a_over_c = 2.0_real64*diffuseness_over_radius
         a_over_c_squared = diffuseness_over_radius*diffuseness_over_radius
         three_a_over_c_squared = 3.0_real64*a_over_c_squared
         a_over_c_cubed = a_over_c_squared*diffuseness_over_radius
         radius_over_diffuseness = half_density_radius/diffuseness
         pi_squared = pi*pi
         half_pi_squared_a_over_c_squared = 0.5_real64*pi_squared*a_over_c_squared
         interior_linear_factor = 1.5_real64 + half_pi_squared_a_over_c_squared
         call es ((-radius_over_diffuseness), origin_inverse_square_sum, origin_inverse_cube_sum)
         six_a_over_c_cubed = 6.0_real64*a_over_c_cubed
         origin_series_correction = -six_a_over_c_cubed*origin_inverse_cube_sum
         fermi_normalization = 1.0_real64 + a_over_c_squared*pi_squared + origin_series_correction
         normalized_nuclear_charge = nuclear_charge/fermi_normalization
         found_first_exterior_point = .false.
         do radial_index = 1, radial_point_count
            radius = radial_grid(radial_index)
            radius_from_half_density = radius - half_density_radius
            scaled_radius_offset = radius_from_half_density/diffuseness
            radius_over_half_density = radius/half_density_radius
            if (radius_over_half_density <= 1.0_real64) then
               call es (scaled_radius_offset, radius_inverse_square_sum, radius_inverse_cube_sum)
               negative_radius_times_potential(radial_index) = normalized_nuclear_charge*(origin_series_correction + &
                  & six_a_over_c_cubed*radius_inverse_cube_sum + radius_over_half_density*(interior_linear_factor - &
                  & three_a_over_c_squared*radius_inverse_square_sum&
                   - 0.5_real64*radius_over_half_density*radius_over_half_density))
            else
               if (.not.found_first_exterior_point) then
                  first_grid_index_beyond_fermi_c = radial_index
                  found_first_exterior_point = .true.
               end if
               call es ((-scaled_radius_offset), radius_inverse_square_sum, radius_inverse_cube_sum)
               negative_radius_times_potential(radial_index) = nuclear_charge*(1.0_real64 + &
                  & three_a_over_c_squared*(radius_over_half_density*radius_inverse_square_sum + &
                  & twice_a_over_c*radius_inverse_cube_sum)/fermi_normalization)
            end if
         end do
      end if
      if (ldbpr(2)) then
         write (99, '(/,"From SUBROUTINE NUCPOT:"/,3(" -------- r -------- ----- -r*V(r) -----"))')
         points_per_column = radial_point_count/3
         if (3*points_per_column == radial_point_count) then
            table_row_count = points_per_column
         else
            table_row_count = points_per_column + 1
         end if
         do table_row = 1, table_row_count
            first_column_index = table_row
            second_column_index = first_column_index + table_row_count
            third_column_index = second_column_index + table_row_count
            if (third_column_index <= radial_point_count) then
               write (99, '(1p,6(1x,1d19.12))') &
                  radial_grid(first_column_index), negative_radius_times_potential(first_column_index), &
                  & radial_grid(second_column_index), &
                  & negative_radius_times_potential(second_column_index), radial_grid(third_column_index), &
                     & negative_radius_times_potential(third_column_index&
                  )
            else if (second_column_index <= radial_point_count) then
               write (99, '(1p,6(1x,1d19.12))') &
                  radial_grid(first_column_index), negative_radius_times_potential(first_column_index), &
                  & radial_grid(second_column_index), &
                  & negative_radius_times_potential(second_column_index)
            else
               write (99, '(1p,6(1x,1d19.12))') &
                  radial_grid(first_column_index), negative_radius_times_potential(first_column_index)
            end if
         end do
      end if
      return
      return
      end subroutine nucpot
end module nuclear_coulomb_potential
