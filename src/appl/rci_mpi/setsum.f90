!***********************************************************************

      subroutine setsum(name)
!
!   Open the  .csum  file on stream 24.
!   Xinghong He                                          10 Jun 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, ierr
!-----------------------------------------------------------------------

      k = index(name,' ')

      call openfl (24, name(1:k-1)//'.csum', 'FORMATTED', 'UNKNOWN', ierr)

      return
      end subroutine setsum
