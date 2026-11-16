!*******************************************************************
!                                                                  *
      integer function jthn(k,n,i)
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vanderbilt University,  Nashville               October  1996  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none (type, external)
      integer, intent(in) :: i, n, k
!-----------------------------------------------
      if(n == 1) then
        jthn=mod(k,i)
      elseif(n == 2) then
        jthn=mod(k/i,i)
      elseif(n == 3) then
        jthn=mod(k/(i*i),i)
      elseif(n == 4) then
        jthn=mod(k/(i*i*i),i)
      elseif(n == 5) then
        jthn=mod(k/(i*i*i*i),i)
      else
        write(6,'(A)') ' ERROR IN JTHN '
        stop
      endif
      return
      end function jthn
