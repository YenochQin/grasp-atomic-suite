module jthn_i
    implicit none (type, external)
    private
    public :: jthn
contains
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
      else if(n == 2) then
        jthn=mod(k/i,i)
      else if(n == 3) then
        jthn=mod(k/(i*i),i)
      else if(n == 4) then
        jthn=mod(k/(i*i*i),i)
      else if(n == 5) then
        jthn=mod(k/(i*i*i*i),i)
      else
        write(6,"(A)") " ERROR IN JTHN "
        stop
      end if
      return
      end function jthn
end module jthn_i
