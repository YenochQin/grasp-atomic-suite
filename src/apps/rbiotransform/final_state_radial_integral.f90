module final_state_radial_integral_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: final_state_radial_integral
contains

    !   The value of RINT is an approximation to:
    !
    !              k
    !         I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)
    !                      i     j       i     j
    !
    !   where   I ( G(r) ; range )  denotes  the  integral  of G(r) over
    !   range.
    !
    ! Written by Farid A Parpia, at Oxford. Last updated: 05 Oct 1992
    real(real64) function final_state_radial_integral (first_orbital, second_orbital, k)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_grid_state, only: radial_grid, radial_grid_derivative
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: final_orbital_last_grid_index, final_radial_large_components, final_radial_small_components
        use radial_quadrature, only: integrate_radial_grid
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_orbital
        integer, intent(in) :: second_orbital
        integer, intent(in) :: k

        ! Local declarations
        integer :: l
        real(real64) :: result
        ! Limit the integral to the supported radial grid points.
        integration_last_index = min(final_orbital_last_grid_index(first_orbital),final_orbital_last_grid_index(second_orbital))
        ! Tabulate the weighted integrand for radial quadrature. The first grid value is unused.
        radial_workspace_a(1) = 0.0_real64
        do l = 2, integration_last_index
            radial_workspace_a(l) = radial_grid(l)**k*(final_radial_large_components(l,first_orbital)* &
                final_radial_large_components(l,second_orbital) + final_radial_small_components(l, &
                first_orbital)*final_radial_small_components(l, &
                & second_orbital))*radial_grid_derivative(l)
        end do
        ! Integrate the tabulated radial integrand.
        call integrate_radial_grid (result)
        final_state_radial_integral = result
        return
    end function final_state_radial_integral
end module final_state_radial_integral_interface
