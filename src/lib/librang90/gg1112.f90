!*******************************************************************
!                                                                  *
      subroutine gg1112(ik1,ik2,bk1,bk2,id1,id2,bd1,bd2,k1,  &
                        qm1,qm2,qm3,qm4,ww)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 15  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!     ELEMENT:                                                     *
!                                                                  *
!       N1       (k1) (j1)(j2)  N1'       N2     (j2)   N2'     +- *
!     (j Q J ::[W(11)*A(1)]  ::j  Q'J')*(j Q J ::A(2)::j Q'J')  -+ *
!       1 1 1                   1  1 1    2 2 2         2 2 2   ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ,WAP1                         *
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
      use cons_c,          only: zero, tenth, half, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use c0t5s_i
      use rmeajj_i
      use wap1_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik1(7), ik2(7), id1(7), id2(7)
      real(real64), intent(in)               :: qm1, qm2, qm3, qm4
      real(real64), intent(in) :: bk1(3), bk2(3), bd1(3), bd2(3)
      real(real64), intent(out)              :: ww
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqmm1, iqmm2, iqmm3, iqmm4, iqmm23, kk1
      real(real64) :: a1, s, bk, wa
!-----------------------------------------------
      ww=zero
      if(ik1(3) > 9) then
        if(ik1(4) > 2) call mes(34)
        if(id1(4) > 2) call mes(34)
      endif
      if(ik2(3) > 9) then
        if(ik2(4) > 2) call mes(34)
        if(id2(4) > 2) call mes(34)
      endif
      iqmm1=qm1+qm1+tenth*qm1
      iqmm2=qm2+qm2+tenth*qm2
      iqmm3=qm3+qm3+tenth*qm3
      iqmm23=iqmm1+iqmm2+iqmm3
      if(ik1(4) /= (id1(4)+iqmm23))return
      iqmm4=qm4+qm4+tenth*qm4
      if(ik2(4) /= (id2(4)+iqmm4))return
      kk1=k1*2
      call c0t5s(bd2(1),bd2(3),qm4,bk2(1),bk2(3),a1)
      if(abs(a1) < eps)return
!GG      CALL SJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S)
      call rmeajj(ik2(3),ik2(1),ik2(7),ik2(6),id2(1),id2(7),id2(6),s)
      if(abs(s) < eps)return
      bk=half*dble(ik2(3))
      call wap1(ik1,bk1,id1,bd1,k1,bk,qm1,qm2,qm3,wa)
      if(abs(wa) < eps)return
      ww=-a1*wa*s/sqrt(dble(ik2(7)+1))
      return
      end subroutine gg1112
