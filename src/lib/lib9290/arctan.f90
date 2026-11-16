!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function arctan (arg1, arg2)
!-----------------------------------------------
!                                                                      *
!                   -1                                                 *
!       ARCTAN = tan   (ARG1/ARG2),      0 .LE. ARCTAN .LT. 2*\pi      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
!   Editted by C. Froese Fischer after translation                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:46:34   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
!
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), intent(in) :: arg1
      real(double), intent(in) :: arg2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical, save :: first = .true., intrin = .true.
!-----------------------------------------------
!
!
!   Determine whether the FORTRAN intrinsic function ATAN2 always
!   returns a positive value
!
      if (first) then
         arctan = atan2(-1.0d00,-1.0d00)
         if (arctan > 0.0d00) then
            intrin = .true.
         else
            intrin = .false.
         endif
         first = .false.
      endif
!
!   Use the intrinsic function if it passes the above test; otherwise
!   add 2*PI to the negative values returned by the intrinsic function
!
      if (intrin) then
         arctan = atan2(arg1,arg2)
      else
         if (arg1 >= 0.0d00) then
            arctan = atan2(arg1,arg2)
         else
            arctan = pi + pi + atan2(arg1,arg2)
         endif
      endif
!
      return
      end function arctan
