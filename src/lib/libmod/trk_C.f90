!*******************************************************************
!                                                                  *
      module trk_c
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 22  --------------   *
!                                                                  *
!   This module is need for librang90.                             *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
!
      implicit none (type, external)
      private :: real64
      integer      :: id1(7),id2(7),ik1(7),ik2(7),id3(7),id4(7),ik3(7),ik4(7)
      real(real64) :: bd1(3),bd2(3),bk1(3),bk2(3),bd3(3),bd4(3),bk3(3),bk4(3)
      end module trk_c
