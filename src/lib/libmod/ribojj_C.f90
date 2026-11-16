!*******************************************************************
!                                                                  *
      module ribojj_c
!                                                                  *
!   This module is need for librang90.                             *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
      implicit none
      integer, dimension(63) :: imptjj, imgtjj, impnjj, imgnjj
      data impnjj/2,1,4,2*3,3*9,3*6,6*18,8*12,20*46,18*26/
      data imgnjj/2,1,5,2*3,3*11,3*8,6*25,8*17,20*63,18*45/
      data imptjj/1,2,3,2*4,3*6,3*9,6*12,8*18,20*26,18*46/
      data imgtjj/1,2,3,2*5,3*8,3*11,6*17,8*25,20*45,18*63/
      end module ribojj_c
