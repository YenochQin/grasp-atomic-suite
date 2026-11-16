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
      use vast_kind_param, only:  double
!
      implicit none
      integer, dimension(7)      :: id1,id2,ik1,ik2,id3,id4,ik3,ik4
      real(double), dimension(3) :: bd1,bd2,bk1,bk2,bd3,bd4,bk3,bk4
!$OMP THREADPRIVATE(ID1, ID2, IK1, IK2, ID3, ID4, IK3, IK4, BD1, BD2, BK1, BK2, BD3, BD4, BK3, BK4)
      end module trk_c
