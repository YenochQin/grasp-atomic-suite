!***********************************************************************
!                                                                      *
      subroutine setmix(name, idblk)
!                                                                      *
!   Opens the  .mix  file on stream 25; writes a header to this file;  *
!   calls LODMIX to interactively determine the eigenpairs required.   *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                                       *
!               [RCI92]: LODMIX.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Modified by Xinghong He               Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      use lodmix_i, only: lodmix
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: name
      character(len=8) :: idblk(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=11), parameter :: form = "UNFORMATTED"
      character(len=7), parameter :: status = "UNKNOWN"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, ierr
!-----------------------------------------------

      k = index(name," ")
      call openfl (25, name(1:k-1)//".cm", form, status, ierr)
      if (ierr /= 0) stop "setmix: Error when opening .cm file"

      write (25) "G92MIX"

      call lodmix (idblk)

      return
      end subroutine setmix
