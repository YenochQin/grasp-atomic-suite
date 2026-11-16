module arctan_i
    implicit none (type, external)
    private
    public :: arctan
contains
!***********************************************************************
!                                                                      *
      real(real64) function arctan (arg1, arg2)
!-----------------------------------------------
!                                                                      *
!                   -1                                                 *
!       ARCTAN = tan   (ARG1/ARG2),      0 .LE. ARCTAN .LT. 2*\pi      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
!   Editted by C. Froese Fischer after translation                     *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
!
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in) :: arg1
      real(real64), intent(in) :: arg2
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
         arctan = atan2(-1.0_real64,-1.0_real64)
         if (arctan > 0.0_real64) then
            intrin = .true.
         else
            intrin = .false.
         end if
         first = .false.
      end if
!
!   Use the intrinsic function if it passes the above test; otherwise
!   add 2*PI to the negative values returned by the intrinsic function
!
      if (intrin) then
         arctan = atan2(arg1,arg2)
      else
         if (arg1 >= 0.0_real64) then
            arctan = atan2(arg1,arg2)
         else
            arctan = pi + pi + atan2(arg1,arg2)
         end if
      end if
!
      return
      end function arctan
end module arctan_i
