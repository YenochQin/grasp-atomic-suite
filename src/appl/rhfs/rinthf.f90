!***********************************************************************
!                                                                      *
      real(real64) function rinthf (i, j, k)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i
      integer, intent(in) :: j
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
      mtp = min(orbital_last_grid_index(i),orbital_last_grid_index(j))
!
!   Value at first tabulation point is arbitrary
!
      ta(1) = 0.0_real64
      do l = 2, mtp
         ta(l) = radial_grid(l)**k*(radial_large_components(l,i)*radial_small_components(l,j) + &
            & radial_small_components(l,i)*radial_large_components(l,j))*radial_grid_derivative(l)
      end do
!
!   Perform integration
!
      call quad (result)
      rinthf = result
!
      return
      end function rinthf
