!***********************************************************************
!                                                                      *
      real(real64) function rinthf (first_orbital, second_orbital, k)
!                                                                      *
!   The value of RINTHF is an approximation to:                        *
!                                                                      *
!              K                                                       *
!         I ( r  * ( P (r)*Q (r) + Q (r)*P (r) ; 0 to infinity)        *
!                     I     J       I     J                            *
!                                                                      *
!   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
!   Range. This is a modification of RINT for the Type B operator.     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                                                      *
!   Written by Per O. Jonsson             Last revision: 24 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use radial_quadrature, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l
      real(real64) :: result
!-----------------------------------------------
!
!   Tabulate integrand as required for SUBROUTINE QUAD
!
      mtp = min(orbital_last_grid_index(first_orbital),orbital_last_grid_index(second_orbital))
!
!   Value at first tabulation point is arbitrary
!
      ta(1) = 0.0_real64
      do l = 2, mtp
         ta(l) = radial_grid(l)**k*(radial_large_components(l,first_orbital)*radial_small_components(l,second_orbital) + &
            & radial_small_components(l,first_orbital)*radial_large_components(l,second_orbital))*radial_grid_derivative(l)
      end do
!
!   Perform integration
!
      call quad (result)
      rinthf = result
!
      return
      end function rinthf
