!***********************************************************************
!                                                                      *
      logical function ldigit (cst)
!                                                                      *
!   .TRUE.  if  CST  is the ASCII representation of a decimal digit;   *
!   .FALSE. otherwise.                                                 *
!                                                                      *
!   Written by Farid A. Parpia             Last revised: 16 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:52   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: cst
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      character, dimension(0:9) :: cdgt
!-----------------------------------------------
!
!
!
      data cdgt/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
!
      do i = 0, 9
         if (cst /= cdgt(i)) cycle
         ldigit = .true.
         go to 2
      end do
      ldigit = .false.
!
    2 continue
      return
      end function ldigit
