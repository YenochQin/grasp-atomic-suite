!***********************************************************************
!                                                                      *
      subroutine fname(name)
!                                                                      *
!   Determines the name of the initial and final states                *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:35:54   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j
!-----------------------------------------------
!
!
!   Obtain the names of the initial and final state files
!
    1 continue
      write (6, *) ' Name of the Initial state'
      read (*, '(A)') name(1)
      write (6, *) ' Name of the Final state'
      read (*, '(A)') name(2)
!
      j = index(name(1),' ')
      if (j == 1) then
         write (6, *) ' Names may not start with blanks'
         go to 1
      endif
!
      j = index(name(2),' ')
      if (j == 1) then
         write (6, *) ' Names may not start with blanks'
         go to 1
      endif

      return
      end subroutine fname
