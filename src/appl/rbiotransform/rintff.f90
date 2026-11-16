!***********************************************************************
!                                                                      *
      real(real64) function rintff (i, j, k)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
!   Determine the maximum tabulation point for the integrand
!
      mtp = min(mfff(i),mfff(j))
!
!   Tabulate the integrand as required for SUBROUTINE QUAD; the
!   value at the first tabulation point is arbitrary
!
      ta(1) = 0.0_real64
      do l = 2, mtp
         ta(l) = radial_grid(l)**k*(pfff(l,i)*pfff(l,j) + qfff(l,i)*qfff(l,j))*radial_grid_derivative(l)
      end do
!
!   Perform the quadrature
!
      call quad (result)
      rintff = result
!
      return
!
      end function rintff
