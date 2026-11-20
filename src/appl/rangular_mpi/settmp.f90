!***********************************************************************
!                                                                      *
      subroutine settmp(nb, kmax, filehead)
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
      integer , intent(in) :: nb
      integer , intent(in) :: kmax
      character , intent(in) :: filehead*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, lck, ierr, i, lng
      character :: ck*2
!-----------------------------------------------
!
! All files  filehead.XX  are UNFORMATTED;
!
      lng = len_trim(filehead)
      do k = 30, 32 + kmax
         call convrt (k, ck, lck)
         call openfl (k, filehead(1:lng)//"."//ck(1:2), "UNFORMATTED", &
            "UNKNOWN", ierr)
         if (ierr == 0) cycle
         do i = 30, k
            close(i)
         end do
         write (stderr_unit, *) "Error when opening the tmp files"
         stop
      end do

      do k = 30, 32 + kmax
         write (k) "MCP", nb
      end do

      return
      end subroutine settmp
