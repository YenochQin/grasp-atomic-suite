module gg1233_i
    implicit none (type, external)
    private
    public :: gg1233
contains
!*******************************************************************
!                                                                  *
      subroutine gg1233(ik1,ik2,ik3,bk1,bk2,bk3,id1,id2,id3,bd1,  &
      bd2,bd3,k1,qm1,qm2,qm3,qm4,ww)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 18  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                      N1     (j1)  N1'      N2     (j2)  N2'      *
!     ELEMENTS:      (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  *
!                      1 1 1        1 1 1    2 2 2        2 2 2    *
!                                                                  *
!        N3     (k1)   N3'                                      +- *
!     *(j Q J::W(33)::j Q'J')                                   -+ *
!        3 3 3         3 3 3                                    ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ,WJ1                          *
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
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use c0t5s_i
      use rmeajj_i
      use wj1_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik1(7),ik2(7),ik3(7),id1(7),id2(7),id3(7)
      real(real64), intent(in)               :: qm1,qm2,qm3,qm4
      real(real64), intent(in) :: bk1(3),bk2(3),bk3(3),bd1(3),bd2(3),bd3(3)
      real(real64), intent(out)              :: ww
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: iqmm1,iqmm2,iqmm3,iqmm4,iqmm34
      real(real64) :: a1, c, c1, s, w
!-----------------------------------------------
      ww=zero
      if(ik1(3) > 9) then
        if(ik1(4) > 2) call mes(32)
        if(id1(4) > 2) call mes(32)
      endif
      if(ik2(3) > 9) then
        if(ik2(4) > 2) call mes(32)
        if(id2(4) > 2) call mes(32)
      endif
      if(ik3(3) > 9) then
        if(ik3(4) > 2) call mes(32)
        if(id3(4) > 2) call mes(32)
      endif
      iqmm1=qm1+qm1+tenth*qm1
      if(ik1(4) /= (id1(4)+iqmm1))return
      iqmm2=qm2+qm2+tenth*qm2
      if(ik2(4) /= (id2(4)+iqmm2))return
      iqmm3=qm3+qm3+tenth*qm3
      iqmm4=qm4+qm4+tenth*qm4
      iqmm34=iqmm3+iqmm4
      if(ik3(4) /= (id3(4)+iqmm34))return
      call c0t5s(bd1(1),bd1(3),qm1,bk1(1),bk1(3),a1)
      if(dabs(a1) < eps)return
      call c0t5s(bd2(1),bd2(3),qm2,bk2(1),bk2(3),c1)
      if(dabs(c1) < eps)return
      a1=a1*c1
      call rmeajj(ik1(3),ik1(1),ik1(7),ik1(6),id1(1),id1(7),id1(6),s)
      if(dabs(s) < eps)return
      call rmeajj(ik2(3),ik2(1),ik2(7),ik2(6),id2(1),id2(7),id2(6),c)
      if(dabs(c) < eps)return
      call wj1(ik3,bk3,id3,bd3,k1,qm3,qm4,w)
      if(dabs(w) < eps)return
      ww=a1*w*s*c/dsqrt(dble((ik1(7)+1)*(ik2(7)+1)))
      return
      end subroutine gg1233
end module gg1233_i
