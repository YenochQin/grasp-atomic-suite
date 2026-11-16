!***********************************************************************
!                                                                      *
      subroutine dampck(ipr, j, ed1, ed2)
!                                                                      *
!   This subroutine determines the damping factor appropriate to the   *
!   present  orbital. The algorithm is taken from C Froese Fischer's   *
!   program MCHF.                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use damp_c
      use orb_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: ipr
      integer, intent(in) :: j
      real(double), intent(inout) :: ed1
      real(double), intent(inout) :: ed2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: adaptv
!-----------------------------------------------
!
!   The damping is adaptive (i.e., can be modified by this SUBROUTINE)
!   if and only if ODAMP(J) .GE. 0.0
!
      adaptv = odamp(j) >= 0.0d00
!
      ed2 = (ed2-e(j))/ed2
      if (adaptv) then
         if (ed1*ed2 < -1.0d-4) then
!           We have a significant oscillation, damp
            odamp(j) = 0.10d00 + 0.90d00*odamp(j)
!           Print *, 'Increasing ODAMP', odamp(j)
         else
            odamp(j) = 0.50d00*odamp(j)
         endif
      endif
!
!     Save the relative difference for the next update of this orbital
      ped(j) = ed2
!
      return
      end subroutine dampck
