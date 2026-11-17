module gg1122_i
    implicit none (type, external)
    private
    public :: gg1122
contains
!*******************************************************************
!                                                                  *
      subroutine gg1122(k1,k2,qm1,qm2,qm3,qm4,aa)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 16  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!     ELEMENTS:                                                    *
!                                                                  *
!       N1       (k1)   N1'         N2       (k2)    N2'        +- *
!     (j  Q J ::W(11)::j  Q'J') * (j  Q J ::W(22):: j  Q'J')    -+ *
!       1  1 1          1  1 1      2  2 2           2  2 2     ++ *
!                                                               -- *
!     SUBROUTINE CALLED: WJ1                                       *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, tenth, eps
      use trk_c, only: bd1, bd2, bd3, bd4, bk1, bk2, bk3, bk4, double, id1, id2, id3, id4, ik1, ik2, ik3, ik4
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i, only: mes
      use wj1_i, only: wj1
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1, k2
      real(real64), intent(in)               :: qm1, qm2, qm3, qm4
      real(real64), intent(out)              :: aa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqmm1, iqmm2, iqmm3, iqmm4, iqmm12, iqmm34
      real(real64) :: a1, w
!-----------------------------------------------
      aa=zero
      if(ik1(3) > 9) then
        if(ik1(4) > 2) call mes(35)
        if(id1(4) > 2) call mes(35)
      end if
      if(ik2(3) > 9) then
        if(ik2(4) > 2) call mes(35)
        if(id2(4) > 2) call mes(35)
      end if
      iqmm1=qm1+qm1+tenth*qm1
      iqmm2=qm2+qm2+tenth*qm2
      iqmm12=iqmm1+iqmm2
      if(ik1(4) /= (id1(4)+iqmm12))return
      iqmm3=qm3+qm3+tenth*qm3
      iqmm4=qm4+qm4+tenth*qm4
      iqmm34=iqmm3+iqmm4
      if(ik2(4) /= (id2(4)+iqmm34))return
      call wj1(ik1,bk1,id1,bd1,k1,qm1,qm2,a1)
      if(dabs(a1) < eps)return
      call wj1(ik2,bk2,id2,bd2,k2,qm3,qm4,w)
      if(dabs(w) < eps)return
      aa=a1*w
      return
      end subroutine gg1122
end module gg1122_i
