!***********************************************************************
!                                                                      *
      subroutine settmp(nb, kmax, filehead)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      use openfl_i
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use iounit_c
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
         call openfl (k, filehead(1:lng)//'.'//ck(1:2), 'UNFORMATTED', &
            'UNKNOWN', ierr)
         if (ierr == 0) cycle
         do i = 30, k
            close(i)
         end do
         write (istde, *) 'Error when opening the tmp files'
         stop
      end do

      do k = 30, 32 + kmax
         write (k) 'MCP', nb
      end do

      return
      end subroutine settmp
