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
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nci
      character, intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: k, ierr
      character :: filnam*256, form*11, status*3
!-----------------------------------------------
!
!   File  rdensity.sum  is FORMATTED
!
      k = index(name,' ')
      if (nci == 0) then
         filnam = name(1:k-1)//'.ci'
      else
         filnam = name(1:k-1)//'.i'
      endif
      form = 'FORMATTED'
      status = 'NEW'
!
      call openfl (24, filnam, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) 'Error when opening', filnam
         stop
      endif
!
      return
      end subroutine setsum
