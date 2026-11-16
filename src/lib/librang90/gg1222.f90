!*******************************************************************
!                                                                  *
      subroutine gg1222(ik1,ik2,bk1,bk2,id1,id2,bd1,bd2,k1, &
                        qm1,qm2,qm3,qm4,ww)
!                                                                  *
! ----------------  SECTION METWO    SUBPROGRAM 17  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                                          N1      (j1)   N1'      *
!     ELEMENTS:                          (j  Q J ::A(1)::j Q'J')*  *
!                                          1  1 1         1 1 1    *
!        N2        (j2)    (k1) (j1)  N2'                       +- *
!     *(j  Q J ::[ A(2) * W(22) ]  ::j  Q'J')                   -+ *
!        2  2 2                       2  2 2                    ++ *
!                                                               -- *
!     SUBROUTINE CALLED: C0T5S,RMEAJJ,AWP1                         *
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
      use awp1_i
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
      integer      :: iqmm1,iqmm2,iqmm3,iqmm4,iqmm34,kk1
      real(real64) :: a1, s, bk, aw
!-----------------------------------------------
      ww=zero
      if(ik1(3) > 9) then
        if(ik1(4) > 2) call mes(33)
        if(id1(4) > 2) call mes(33)
      endif
      if(ik2(3) > 9) then
        if(ik2(4) > 2) call mes(33)
        if(id2(4) > 2) call mes(33)
      endif
      iqmm1=qm1+qm1+tenth*qm1
      if(ik1(4) /= (id1(4)+iqmm1))return
      iqmm2=qm2+qm2+tenth*qm2
      iqmm3=qm3+qm3+tenth*qm3
      iqmm4=qm4+qm4+tenth*qm4
      iqmm34=iqmm2+iqmm3+iqmm4
      if(ik2(4) /= (id2(4)+iqmm34))return
      kk1=k1*2
      call c0t5s(bd1(1),bd1(3),qm1,bk1(1),bk1(3),a1)
      if(dabs(a1) < eps)return
!GG      CALL SJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      call rmeajj(ik1(3),ik1(1),ik1(7),ik1(6),id1(1),id1(7),id1(6),s)
      if(dabs(s) < eps)return
      bk=half*dble(ik1(3))
      call awp1(ik2,bk2,id2,bd2,k1,bk,qm2,qm3,qm4,aw)
      if(dabs(aw) < eps)return
      ww=-a1*aw*s/dsqrt(dble(ik1(7)+1))
      return
      end subroutine gg1222
