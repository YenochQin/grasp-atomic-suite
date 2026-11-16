!***********************************************************************
!                                                                      *
      subroutine fname(name)
!                                                                      *
!   Determines the name of the initial and final states                *
!   In addition this subroutine determines which J symmetries          *
!   that are to be transformed                                         *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:19:37   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   C O M M O N   B l o c k s
!-----------------------------------------------
      use jqjc_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i
      logical :: yes
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
! Per april 2007
! Check if the initial and final states are identical.

      if (trim(name(1)) == trim(name(2))) then
         print *
         print *, ' Initial and final states are identical and there is'
         print *, ' no need for the biorthogonal transformation. Just  '
         print *, ' copy name.w to name.bw and name.(c)m to name.(c)bm '
         print *, ' and run rtransition.                               '
         print *


         print *, ' Do you want to continue anyway ? '
         yes = getyn ()
         if (yes.eqv..false.) stop
      end if

! end Per 2007
      write (6, *) ' Transformation of all J symmetries?'
      yes = getyn()
      if (yes) then
         ntrans = 0
      else
         ntrans = 1
         write (6, *) ' Number of initial state J symmetries to be transformed'
         read (*, *) jqj1
         write (6, *) ' Give the J symmetries in the form 2*J'
         read (*, *) (itjqj1(i),i=1,jqj1)
         itjqj1(:jqj1) = itjqj1(:jqj1) + 1
         write (6, *) ' Number of final state J symmetries to be transformed'
         read (*, *) jqj2
         write (6, *) ' Give the J symmetries in the form 2*J'
         read (*, *) (itjqj2(i),i=1,jqj2)
         itjqj2(:jqj2) = itjqj2(:jqj2) + 1
      endif

      return
      end subroutine fname
