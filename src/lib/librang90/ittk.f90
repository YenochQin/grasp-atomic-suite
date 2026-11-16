module ittk_i
    implicit none (type, external)
    private
    public :: ittk
contains
!*******************************************************************
!                                                                  *
      integer function ittk (i, j, k)
!                                                                  *
!     CHESKED TRIANGULAR CONDITIONS FOR   I/2, J/2, K/2.           *
!     I+J>=K, I+K>=J, J+K>=I,                                      *
!     I/2+J/2+K/2 - WHOLE NUMBER                                   *
!     ITTK=1 -   IF NOT SATISFY                                    *
!     ITTK=0 -   IN OVER CASES                                     *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i, j, k
!-----------------------------------------------
      ittk = 0
      if (iabs(i - j) > k) return
      if (i + j < k) return
      if (mod(i + j + k,2) /= 0) return
      ittk = 1
      return
      end function ittk
end module ittk_i
