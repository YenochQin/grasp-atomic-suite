!***********************************************************************
!                                                                      *
      subroutine dampck(ipr, j, ed1, ed2)
!                                                                      *
!   This subroutine determines the damping factor appropriate to the   *
!   present  orbital. The algorithm is taken from C Froese Fischer's   *
!   program MCHF.                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!   Modified by C. Froese Fischer           Last update: 07 Apr 2009   *
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
      if (ipr /= j) then
         if (adaptv) odamp(j) = 0.75d00*odamp(j)
      else
         ed2 = ed2-e(j)
         if (adaptv) then
            if (ed1*ed2 > 0.0d00) then
               odamp(j) = 0.75d00*odamp(j)
            else
               odamp(j) = 0.25d00+0.75d00*odamp(j)
            endif
         endif
      endif
!
      if (ipr == j) then
         ed1 = ed2
      else
         ed1 = ed2-e(j)
      endif
!
      ipr = j
      return
      end subroutine dampck
