module getyn_i
    implicit none (type, external)
    private
    public :: getyn
contains
!***********************************************************************
!                                                                      *
      logical function getyn ()
!                                                                      *
!   This  subprogram reads a response on  the default input unit; the  *
!   response must be either 'y' or 'n'. GETYN is .TRUE. if 'y' is en-  *
!   tered and .FALSE. if 'n' is entered.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 27 Aug 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c,     only: stdin_unit, stderr_unit
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character (len = 1) :: rspns
!-----------------------------------------------
    1 continue
      read (stdin_unit, "(A)") rspns
      if (rspns == "y" .or. rspns == "Y") then
         getyn = .true.
      else if (rspns == "n" .or. rspns == "N") then
         getyn = .false.
      else
         write (stderr_unit, *) "Expecting <y><cr> or <n><cr> ..."
         goto 1
      end if
!
      return
      end function getyn
end module getyn_i
