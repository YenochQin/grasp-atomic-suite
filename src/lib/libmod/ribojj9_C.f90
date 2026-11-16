!*******************************************************************
!                                                                  *
      module ribojj9_c
!                                                                  *
!   This module is need for librang90.                             *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
      implicit none
      integer, dimension(6) :: imptjj9, imgtjj9, impnjj9, imgnjj9
      data imptjj9/301,5*302/
      data imgtjj9/301,5*306/
      data impnjj9/302,5*301/
      data imgnjj9/306,5*301/
      end module ribojj9_c
