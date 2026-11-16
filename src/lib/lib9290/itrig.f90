!***********************************************************************
!                                                                      *
      integer function itrig (i1, i2, i3)
!                                                                      *
!   The  triangular delta. Input: Values of 2*J+1; Output: 1, IF J'S   *
!   form a triangle; 0, otherwise.                                     *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:47   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i1
      integer, intent(in) :: i2
      integer, intent(in) :: i3
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i4
!-----------------------------------------------
!
      i4 = i2 - i3
      if (i1>=abs(i4) + 1 .and. i1<=i2+i3-1) then
         itrig = 1
      else
         itrig = 0
      endif
!
      return
      end function itrig
