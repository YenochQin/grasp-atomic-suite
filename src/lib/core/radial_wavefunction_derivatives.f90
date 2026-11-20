! Tabulate step-scaled derivatives of both radial spinor components on the transformed grid.
module radial_wavefunction_derivatives
    implicit none (type, external)
    private
    public :: dpbdt
contains
      subroutine dpbdt(orbital_index)
      use, intrinsic :: iso_fortran_env, only: real64
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_derivative_stencil, only: derivative_weights => a
      use radial_quadrature_workspace, only: scaled_large_derivative_grid => ta, scaled_small_derivative_grid => tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer :: radial_index, k, stencil_row, loc
      real(real64) :: offset_6_weight, offset_5_weight, offset_4_weight, offset_3_weight, offset_2_weight, &
         & offset_1_weight, scaled_large_derivative, scaled_small_derivative, aik

      offset_6_weight=derivative_weights(7,1)
      offset_5_weight=derivative_weights(7,2)
      offset_4_weight=derivative_weights(7,3)
      offset_3_weight=derivative_weights(7,4)
      offset_2_weight=derivative_weights(7,5)
      offset_1_weight=derivative_weights(7,6)
      ! Use one-sided rows at both edges and the antisymmetric centered row in the interior.
      do radial_index = 1,6
         scaled_large_derivative = 0.0_real64
         scaled_small_derivative = 0.0_real64
         scaled_large_derivative = sum(derivative_weights(radial_index   ,:)*radial_large_components(:13,orbital_index))
         scaled_small_derivative = sum(derivative_weights(radial_index   ,:)*radial_small_components(:13,orbital_index))
         scaled_large_derivative_grid(radial_index) = scaled_large_derivative
         scaled_small_derivative_grid(radial_index) = scaled_small_derivative
      end do

      do radial_index = 7, radial_point_count - 6
         scaled_large_derivative_grid(radial_index) = offset_6_weight*(radial_large_components(radial_index - 6, &
            & orbital_index) - radial_large_components(radial_index + 6,orbital_index)) + offset_5_weight* &
            & (radial_large_components(radial_index - 5,orbital_index) - radial_large_components(radial_index + 5, &
               & orbital_index)&
            ) + offset_4_weight*(radial_large_components(radial_index - 4, &
               & orbital_index) - radial_large_components(radial_index + 4,orbital_index)) + offset_3_weight* &
               & (radial_large_components(radial_index - 3,orbital_index) - radial_large_components(radial_index + 3, &
                  & orbital_index))&
             + offset_2_weight*(radial_large_components(radial_index - 2, &
                & orbital_index) - radial_large_components(radial_index + 2,orbital_index)) + offset_1_weight* &
                & (radial_large_components(radial_index - 1,orbital_index) - radial_large_components(radial_index + 1, &
                   & orbital_index))
         scaled_small_derivative_grid(radial_index) = offset_6_weight*(radial_small_components(radial_index - 6, &
            & orbital_index) - radial_small_components(radial_index + 6,orbital_index)) + offset_5_weight* &
            & (radial_small_components(radial_index - 5,orbital_index) - radial_small_components(radial_index + 5, &
               & orbital_index)&
            ) + offset_4_weight*(radial_small_components(radial_index - 4, &
               & orbital_index) - radial_small_components(radial_index + 4,orbital_index)) + offset_3_weight* &
               & (radial_small_components(radial_index - 3,orbital_index) - radial_small_components(radial_index + 3, &
                  & orbital_index))&
             + offset_2_weight*(radial_small_components(radial_index - 2, &
                & orbital_index) - radial_small_components(radial_index + 2,orbital_index)) + offset_1_weight* &
                & (radial_small_components(radial_index - 1,orbital_index) - radial_small_components(radial_index + 1, &
                   & orbital_index))
      end do
      do radial_index = radial_point_count - 5, radial_point_count
         stencil_row = radial_index - radial_point_count + 13
          scaled_large_derivative = sum(derivative_weights(stencil_row, &
             & :)*radial_large_components(radial_point_count-12:radial_point_count,orbital_index))
          scaled_small_derivative = sum(derivative_weights(stencil_row, &
             & :)*radial_small_components(radial_point_count-12:radial_point_count,orbital_index))
         scaled_large_derivative_grid(radial_index) = scaled_large_derivative
         scaled_small_derivative_grid(radial_index) = scaled_small_derivative
      end do
      return
      end subroutine dpbdt
end module radial_wavefunction_derivatives
