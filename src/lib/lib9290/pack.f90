module pack_implementation
    implicit none (type, external)
    private
    public :: pack_impl
contains
!***********************************************************************
!                                                                      *
      subroutine pack_impl(iunpkd, isubsh, ipackd)
!                                                                      *
!   Subshell occupation numbers and all angular momenta 2J+1 are not   *
!   likely to exceed 127 in any application of the GRASP92 suite. It   *
!   is, therefore, inefficient to allocate an entire INTEGER storage   *
!   cell to  any of these quantities --- a single  byte is adequate.   *
!   Up to eight integers of magnitude less than or equal to  127 may   *
!   be stored in one  64-bit cell, four in a 32-bit cell.  This idea   *
!   is implemented in  the present subprogram.  IPACKD is assumed to   *
!   be an INTEGER vector of at least NINT (NW/8) elements for 64-bit   *
!   architectures, and NINT (NW/4) elements for 32-bit architectures.  *
!   ISUBSH  is the subshell sequence number, IUNPKD its value.  LOC1   *
!   is the element number;  LOC2 is one less than the byte number in   *
!   this element.                                                      *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
!   Modified by G. Gaigalas                                 May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: int8
      use iounit_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: iunpkd
      integer, intent(in) :: isubsh
      integer(int8), intent(inout) :: ipackd(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
      if (abs(iunpkd) > 127) then
         write (stderr_unit, *) "PACK: Argument IUNPKD out of range."
         stop
      end if
!
      if (isubsh <= 0) then
         write (stderr_unit, *) "PACK: ISUBSH=", isubsh, " less than 1"
         stop
      end if
!
      ipackd(isubsh) = iunpkd
!
      return
      end subroutine pack_impl
end module pack_implementation

! GNU Fortran 15 confuses intrinsic names with module procedure pointers.
! Re-export the implementation under the established public name.
module pack_i
    use pack_implementation, only: pack => pack_impl
    implicit none (type, external)
    private
    public :: pack
end module pack_i
