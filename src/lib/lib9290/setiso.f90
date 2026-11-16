module setiso_i
    implicit none (type, external)
    private
    public :: setiso
contains
!***********************************************************************
!                                                                      *
      subroutine setiso(fname)
!                                                                      *
!   Open, check, load data from and close the  .iso  file. This file   *
!   is always attached to stream  22 in  RSCF92,  RCI92,  HFS92, and   *
!   OSCL92.                                                            *
!   Filename - fname is moved to the argument list. This subroutine
!   is ready to replace the one in lib/lib92, but need to check
!   other application programs before doing do.
!                                                                      *
!   Call(s) to: [LIB92]: LODISO, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified    by Xinghong He            Last revision:  1 Jun 1998   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      use lodiso_i, only: lodiso
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = *), intent(in) :: fname
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=3), parameter  :: status = "OLD"
      character(len=6), parameter  :: myname = "SETISO"
      character(len=9), parameter  :: form = "FORMATTED"
      character(len=14), parameter :: signature = "Atomic number:"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer           :: lenf, ierr, ios
      logical           :: found
      character(len=14) :: str
!-----------------------------------------------
!
      inquire(file=fname, exist=found)
      if (.not.found) then
         lenf = len_trim(fname)
         write (stderr_unit,*)                                              &
                     myname,"- file: ",fname(1:lenf)," does not exist"
         stop
      end if
!
      call openfl (22, fname, form, status, ierr)
      if (ierr /= 0) then
         write (stderr_unit, *) "Error opening isodata file: ", fname(1:lenf)
         stop
      end if
!
!   Check the first record of the file; if not as expected, try again
!
      read (22, "(A)", iostat=ios) str
      if (ios/=0 .or. str/=signature) then
         write (stderr_unit, *) "Not an ISOtope Data File;"
         close(22)
         stop
      end if
!
!   Load data from the .iso file and then close it.
!
      call lodiso
      close(22)

      return
      end subroutine setiso
end module setiso_i
