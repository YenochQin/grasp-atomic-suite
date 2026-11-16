module setcsh_i
    implicit none (type, external)
    private
    public :: setcsh
contains
!***********************************************************************
!
      subroutine setcsh(nfile, name, ncore)
!
!   Open, check the CSL file and load the load (via lodcsh) data from
!   the header lines. It is designed to replace all kinds of "setcsl"
!   routines within GRASP packages.
!
!   Routines called: lodcsh
!
!   ncore is output parameter
!
!   Written by Xinghone He                Last revision: 23 Dec 1997
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lodcsh_i, only: lodcsh
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer, intent(out) :: ncore
      character (len = *), intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: length, ios
      logical :: found
      character :: record*15
!-----------------------------------------------
!
!     ...Locals

      length = len_trim(name)

      inquire(file=name, exist=found)
      if (.not.found) then
         write (stderr_unit, *) name(1:length), " does not exist"
         stop
      end if

      inquire(unit=nfile, opened=found)
      if (found) then
         write (stderr_unit, *) "Unit ", nfile, " has been used elsewhere"
         stop
      end if

      open(nfile, file=name, status="OLD", iostat=ios, position="asis")
      if (ios /= 0) then
         write (stderr_unit, *) "Error when opening ", name
         stop
      end if
!
!   Check the first record of the file; if not as expected, try again
!
      read (nfile, "(1A15)", iostat=ios) record

      if (ios/=0 .or. record(1:15)/="Core subshells:") then
         write (stderr_unit, *) "Not a Configuration Symmetry List File;"
         close(nfile)
         stop
      end if
!
!   Load data from the  .csl  file
!
      call lodcsh (nfile, ncore)

      return
      end subroutine setcsh
end module setcsh_i
