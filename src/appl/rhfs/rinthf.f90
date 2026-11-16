!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function rinthf (i, j, k)
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
      use vast_kind_param, only:  double
      use grid_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i
      implicit none
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
!   Tabulate integrand as required for SUBROUTINE QUAD
!
      mtp = min(mf(i),mf(j))
!
!   Value at first tabulation point is arbitrary
!
      ta(1) = 0.0d00
      do l = 2, mtp
         ta(l) = r(l)**k*(pf(l,i)*qf(l,j) + qf(l,i)*pf(l,j))*rp(l)
      end do
!
!   Perform integration
!
      call quad (result)
      rinthf = result
!
      return
      end function rinthf
