!***********************************************************************
!                                                                      *
      subroutine settmpgg(nb, k, filehead)
!***********************************************************************
!   Modified by Gediminas Gaigalas:                         Feb 2017   *
!                               1) for new spin-angular integration,   *
!                               2) for sorting in the memory.          *
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      use file_open, only: openfl
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nb
      integer, intent(in) :: k
      character, intent(in) :: filehead*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: lck, ierr, lng
      character :: ck*2
!-----------------------------------------------
!
! All files  filehead.XX  are UNFORMATTED;
!
      lng = len_trim(filehead)
      call convrt (k, ck, lck)
      call openfl (k, filehead(1:lng)//"."//ck(1:2), "UNFORMATTED", &
            "UNKNOWN", ierr)
      if (ierr /= 0) then
         close(k)
         write (stderr_unit, *) "Error when opening the tmp files"
         stop
      end if
      write (k) "MCP", nb
      return
      end subroutine settmpgg
