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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
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
      character :: filnam1*256, filnam2*256, defnam*11, form*11, status*3
!-----------------------------------------------
!
!   File  hfs92.sum  is FORMATTED
!
      k = index(name," ")
      if (nci == 0) then
         filnam1 = name(1:k-1)//".ch"
         filnam2 = name(1:k-1)//".choffd"
      else
         filnam1 = name(1:k-1)//".h"
         filnam2 = name(1:k-1)//".hoffd"
      end if
      form = "FORMATTED"
      status = "NEW"
!
      call openfl (29, filnam1, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) "Error when opening", filnam1
         stop
      end if
!
      call openfl (24, filnam2, form, status, ierr)
      if (ierr /= 0) then
         write (6, *) "Error when opening", filnam2
         stop
      end if
!
      return
      end subroutine setsum
