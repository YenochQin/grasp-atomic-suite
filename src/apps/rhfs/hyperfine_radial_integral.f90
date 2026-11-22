! Compute the radial part of a relativistic hyperfine matrix element.
module hyperfine_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: hyperfine_radial_integral
contains

    !   The value of RINTHF is an approximation to:
    !
    !              K
    !         I ( r  * ( P (r)*Q (r) + Q (r)*P (r) ; 0 to infinity)
    !
    !   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over
    !   Range. This is a modification of RINT for the Type B operator.
    !
    ! Written by Per O. Jonsson. Last revision: 24 Dec 1992
    real(real64) function hyperfine_radial_integral (first_orbital, second_orbital, radial_power)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_grid, radial_grid_derivative
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital
        integer, intent(in) :: radial_power

        ! Local declarations
        integer :: radial_index
        real(real64) :: radial_integral
        ! Tabulate integrand as required for SUBROUTINE QUAD
        integration_last_index = min(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
        ! Value at first tabulation point is arbitrary
        radial_workspace_a(1) = 0.0_real64
        do radial_index = 2, integration_last_index
            radial_workspace_a(radial_index) = radial_grid(radial_index)**radial_power*(radial_large_components(radial_index, &
                & first_orbital)*radial_small_components(radial_index, &
                second_orbital) + &
                & radial_small_components(radial_index,first_orbital)*radial_large_components(radial_index,second_orbital))* &
                & radial_grid_derivative(radial_index)
        end do
        ! Perform integration
        call integrate_radial_grid (radial_integral)
        hyperfine_radial_integral = radial_integral
        return
    end function hyperfine_radial_integral
end module hyperfine_radial_integral_interface
