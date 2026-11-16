!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function rint (i, j, k)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:22   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use grid_c
      use tatb_c, only: ta, mtp
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i
      integer  :: j
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l
      real(double) :: result
!-----------------------------------------------
!
!
!   Determine the maximum tabulation point for the integrand
!
      mtp = min(orbital_last_grid_index(i),orbital_last_grid_index(j))
!
!   Tabulate the integrand as required for SUBROUTINE QUAD; the
!   value at the first tabulation point is arbitrary
!
      ta(1) = 0.0d00
      do l = 2, mtp
         ta(l) = radial_grid(l)**k*(radial_large_components(l,i)*radial_large_components(l,j) + &
            & radial_small_components(l,i)*radial_small_components(l,j))*radial_grid_derivative(l)
      end do
!
!   Perform the quadrature
!
      call quad (result)
      rint = result
!
      return
!
      end function rint
