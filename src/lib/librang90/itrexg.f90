!*******************************************************************
!                                                                  *
      integer function itrexg(i1,i2,i3,i4,k)
!
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: i1, i2, i3, i4
      integer, intent(out) :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j
!-----------------------------------------------
      j=max0(iabs(i1-i2),iabs(i3-i4))
      k=min0(iabs(i1+i2),iabs(i3+i4))-j+1
      itrexg=j
      return
      end function itrexg
