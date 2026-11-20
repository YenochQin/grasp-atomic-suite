!***********************************************************************

      subroutine setsum(name)
!
!   Open the  .csum  file on stream 24.
!   Xinghong He                                          10 Jun 1998
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr
!-----------------------------------------------------------------------


      call openfl (24, trim(name)//".csum", "FORMATTED", "UNKNOWN", ierr)

      return
      end subroutine setsum
