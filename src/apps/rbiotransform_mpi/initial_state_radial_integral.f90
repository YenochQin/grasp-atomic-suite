!***********************************************************************
!                                                                      *
real(real64) function initial_state_radial_integral (first_orbital, second_orbital, k)
    ! Module dependencies
    !                                                                      *
    !   The value of RINT is an approximation to:                          *
    !                                                                      *
    !              k                                                       *
    !         I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)       *
    !                      i     j       i     j                           *
    !                                                                      *
    !   where   I ( G(r) ; range )  denotes  the  integral  of G(r) over   *
    !   range.                                                             *
    !                                                                      *
    !   Call(s) to: [LIB92]: QUAD.                                         *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford   Last updated: 05 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use radial_quadrature, only: quad
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: first_orbital
    integer, intent(in) :: second_orbital
    integer, intent(in) :: k

    ! Local declarations
    integer :: l
    real(real64) :: result
    !-----------------------------------------------
    !
    !   Determine the maximum tabulation point for the integrand
    !
    mtp = min(mfii(first_orbital),mfii(second_orbital))
    !
    !   Tabulate the integrand as required for SUBROUTINE QUAD; the
    !   value at the first tabulation point is arbitrary
    !
    ta(1) = 0.0_real64
    do l = 2, mtp
        ta(l) = radial_grid(l)**k*(pfii(l,first_orbital)*pfii(l,second_orbital) + qfii(l,first_orbital)*qfii(l, &
            & second_orbital))*radial_grid_derivative(l)
    end do
    !
    !   Perform the quadrature
    !
    call quad (result)
    initial_state_radial_integral = result
    !
    return
    !
end function initial_state_radial_integral
