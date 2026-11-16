!*******************************************************************
!                                                                  *
      integer function ixjtik (i, j, k, l, m, n)
!                                                                  *
!     CHESKED TRIANGULAR CONDITIONS FOR 6j COEFFICIENT             *
!                                                                  *
!     | I/2  J/2  K/2 |            IXJTIK=1 - IF NOT SATISFY       *
!     | L/2  M/2  N/2 |            IXJTIK=0 - IN OVER CASES        *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i, j, k, l, m, n
!-----------------------------------------------
      ixjtik = 0
      if (ittk(i,j,k) == 0) return
      if (ittk(i,m,n) == 0) return
      if (ittk(l,j,n) == 0) return
      if (ittk(l,m,k) == 0) return
      ixjtik = 1
      return
      end function ixjtik
