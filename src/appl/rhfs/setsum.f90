!***********************************************************************
!                                                                      *
      subroutine setsum(name, nci)
!                                                                      *
!   Open the  .sum  files on stream 24 and 29                          *
!                                                                      *
!   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!                                                                      *
!   Updated by Per Jonsson                               28 Oct 1999   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nci
      character, intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: k, ierr
      character :: filnam1*256, filnam2*256, defnam*11, form*11, status*3
!-----------------------------------------------
!
!   File  hfs92.sum  is FORMATTED
!
      k = index(name,' ')
      if (nci == 0) then
         filnam1 = name(1:k-1)//'.ch'
         filnam2 = name(1:k-1)//'.choffd'
      else
         filnam1 = name(1:k-1)//'.h'
         filnam2 = name(1:k-1)//'.hoffd'
      endif
      form = 'FORMATTED'
      status = 'NEW'
!
      call openfl (29, filnam1, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) 'Error when opening', filnam1
         stop
      endif
!
      call openfl (24, filnam2, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) 'Error when opening', filnam2
         stop
      endif
!
      return
      end subroutine setsum
