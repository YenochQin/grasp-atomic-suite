! Tabulate step-scaled derivatives of both radial spinor components on the transformed grid.
module radial_wavefunction_derivatives
    implicit none (type, external)

    ! Module interface
    private
    public :: differentiate_radial_wavefunctions
contains
    subroutine differentiate_radial_wavefunctions(orbital_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_point_count
        use radial_derivative_stencil, only: derivative_weights
        use radial_quadrature_workspace, only: radial_workspace_a, radial_workspace_b
        use radial_wavefunctions, only: radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index

        ! Local declarations
        integer :: radial_index, stencil_row
        real(real64) :: offset_6_weight, offset_5_weight, offset_4_weight, offset_3_weight, offset_2_weight, &
            & offset_1_weight, scaled_large_derivative, scaled_small_derivative

        ! The unit-spaced stencil stores h*dP/dt in workspace A and h*dQ/dt in workspace B.
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
            radial_workspace_a(radial_index) = scaled_large_derivative
            radial_workspace_b(radial_index) = scaled_small_derivative
        end do

        do radial_index = 7, radial_point_count - 6
            radial_workspace_a(radial_index) = offset_6_weight*(radial_large_components(radial_index - 6, &
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
            radial_workspace_b(radial_index) = offset_6_weight*(radial_small_components(radial_index - 6, &
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
            radial_workspace_a(radial_index) = scaled_large_derivative
            radial_workspace_b(radial_index) = scaled_small_derivative
        end do
        return
    end subroutine differentiate_radial_wavefunctions
end module radial_wavefunction_derivatives
