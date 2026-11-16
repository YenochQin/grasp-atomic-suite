!***********************************************************************
      subroutine setsum(filnam)
!***********************************************************************
!-----------------------------------------------
!    M O D U L E S
!-----------------------------------------------
     use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: filnam*(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=9), parameter :: form = 'FORMATTED'
      character(len=3), parameter :: status = 'NEW'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
!-----------------------------------------------

      call openfl (24, filnam, form, status, ierr)
      if (ierr /= 0) then
         write (stderr_unit, *) 'Error when opening ', filnam
         stop
      endif

      return
      end subroutine setsum
