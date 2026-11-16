module gg1234_i
    implicit none (type, external)
    private
    public :: gg1234
contains
!*******************************************************************
!                                                                  *
      subroutine gg1234(ik1,ik2,ik3,ik4,bk1,bk2,bk3,bk4,id1,id2,  &
      id3,id4,bd1,bd2,bd3,bd4,qm1,qm2,qm3,qm4,ww)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 19  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                      N1     (j1)  N1'      N2    (j2)   N2'      *
!     ELEMENTS:      (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  *
!                      1 1 1        1 1 1    2 2 2        2 2 2    *
!                                                                  *
!        N3     (j3)  N3'      N4     (j4)  N4'                 +- *
!     *(j Q J::A(3)::j Q'J')*(j Q J::A(4)::j Q'J')              -+ *
!        3 3 3        3 3 3    4 4 4        4 4 4               ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ                              *
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
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in) :: ik1(7),ik2(7),ik3(7),ik4(7),id1(7),id2(7),id3(7),id4(7)
      real(real64), intent(in)               :: qm1,qm2,qm3,qm4
      real(real64), intent(in) :: bk1(3),bk2(3),bk3(3),bk4(3),bd1(3),bd2(3),bd3(3),bd4(3)
      real(real64), intent(out)              :: ww
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqmm1, iqmm2, iqmm3, iqmm4
      real(real64) :: a1, c, c1, s, v, z
!-----------------------------------------------
      ww=zero
      if(ik1(3) > 9) then
        if(ik1(4) > 2) call mes(30)
        if(id1(4) > 2) call mes(30)
      endif
      if(ik2(3) > 9) then
        if(ik2(4) > 2) call mes(30)
        if(id2(4) > 2) call mes(30)
      endif
      if(ik3(3) > 9) then
        if(ik3(4) > 2) call mes(30)
        if(id3(4) > 2) call mes(30)
      endif
      if(ik4(3) > 9) then
        if(ik4(4) > 2) call mes(30)
        if(id4(4) > 2) call mes(30)
      endif
      iqmm1=qm1+qm1+tenth*qm1
      if(ik1(4) /= (id1(4)+iqmm1))return
      iqmm2=qm2+qm2+tenth*qm2
      if(ik2(4) /= (id2(4)+iqmm2))return
      iqmm3=qm3+qm3+tenth*qm3
      if(ik3(4) /= (id3(4)+iqmm3))return
      iqmm4=qm4+qm4+tenth*qm4
      if(ik4(4) /= (id4(4)+iqmm4))return
      call c0t5s(bd1(1),bd1(3),qm1,bk1(1),bk1(3),a1)
      if(dabs(a1) < eps)return
      call c0t5s(bd2(1),bd2(3),qm2,bk2(1),bk2(3),c1)
      if(dabs(c1) < eps)return
      a1=a1*c1
      call c0t5s(bd3(1),bd3(3),qm3,bk3(1),bk3(3),c1)
      if(dabs(c1) < eps)return
      a1=a1*c1
      call c0t5s(bd4(1),bd4(3),qm4,bk4(1),bk4(3),c1)
      if(dabs(c1) < eps)return
      a1=a1*c1
!GG      CALL SJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      call rmeajj(ik1(3),ik1(1),ik1(7),ik1(6),id1(1),id1(7),id1(6),s)
      if(dabs(s) < eps)return
!GG      CALL SJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
      call rmeajj(ik2(3),ik2(1),ik2(7),ik2(6),id2(1),id2(7),id2(6),c)
      if(dabs(c) < eps)return
!GG      CALL SJJ(IK3(3),IK3(1),IK3(7),IK3(6),ID3(1),ID3(7),ID3(6),V)
      call rmeajj(ik3(3),ik3(1),ik3(7),ik3(6),id3(1),id3(7),id3(6),v)
      if(dabs(v) < eps)return
!GG      CALL SJJ(IK4(3),IK4(1),IK4(7),IK4(6),ID4(1),ID4(7),ID4(6),Z)
      call rmeajj(ik4(3),ik4(1),ik4(7),ik4(6),id4(1),id4(7),id4(6),z)
      if(dabs(z) < eps)return
      ww=a1*s*c*v*z/dsqrt(dble((ik1(7)+1)*(ik2(7)+1)*(ik3(7)+1)*(ik4(7)+1)))
      return
      end subroutine gg1234
end module gg1234_i
