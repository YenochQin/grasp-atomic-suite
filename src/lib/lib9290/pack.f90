!***********************************************************************
!                                                                      *
      subroutine pack(iunpkd, isubsh, ipackd)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:08   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  byte
      use iounit_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: iunpkd
      integer, intent(in) :: isubsh
      integer(byte), dimension(*), intent(inout) :: ipackd
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
      if (abs(iunpkd) > 127) then
         write (istde, *) 'PACK: Argument IUNPKD out of range.'
         stop
      endif
!
      if (isubsh <= 0) then
         write (istde, *) 'PACK: ISUBSH=', isubsh, ' less than 1'
         stop
      endif
!
      ipackd(isubsh) = iunpkd
!
      return
      end subroutine pack
