!*******************************************************************
!                                                                  *
      subroutine wj1(ik,bk,id,bd,k2,qm1,qm2,wj)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 24  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                      N       (k2)  N'     +-                     *
!     ELEMENT:       (j  QJ:: W   ::j  QJ)  -+                     *
!                                           ++                     *
!                                           -- S5(1.47),(1.48),    *
!                                                (1.49),(1.50).    *
!                                                                  *
!     SUBROUTINE CALLED: CLE0SM,C1E1SM,RWJJ                        *
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
      use vast_kind_param, only: double
      use cons_c,          only: zero, tenth, one, two, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use w1jjg_i
      use cle0sm_i
      use c1e1sm_i
      use rwjj_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k2
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in)               :: qm1, qm2
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(out)              :: wj
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: k1, iq, iqm2
      real(double) :: a, qq, w, wk1
!-----------------------------------------------
      wj=zero
      if(id(3) == 9) then
        if(max0(ik(4),id(4)) < 3) then
          if(ik(1) < 300) call mes(56)
          if(id(1) < 300) call mes(56)
          iqm2=qm2+qm2+qm2*eps
          if((id(4)+iqm2) > 2) call mes(2)
          call w1jjg(k2,qm1,qm2,ik,bk,id,bd,wj)
          return
        else
          print*, "ERROR in  WJ1"
          stop
        endif
      elseif(id(3) > 9) then
        iqm2=qm2+qm2+qm2*eps
        if((id(4)+iqm2) > 2) call mes(2)
        call w1jjg(k2,qm1,qm2,ik,bk,id,bd,wj)
        return
      endif
      qq=qm1+qm2
      if(dabs(qq) >= eps) then
        if(((k2/2)*2) /= k2)return
        iq=qq+qq*tenth
        if((ik(4)-id(4)-2*iq) /= 0)return
        call c1e1sm(bd(1),bd(3),qq,bk(1),bk(3),a)
        if(dabs(a) < eps)return
        call rwjj(ik(3),ik(1),id(1),1,k2,w)
        wj=a*w/dsqrt(two*bk(1)+one)
      else if(ik(4) /= id(4))then
        return
      else if(k2 == 0) then
        if(id(1) /= ik(1))return
        if(qm1 < eps) then
          a=dble(id(3)+1-id(4))
        else
          a=-dble(id(4))
        end if
        wj=a*dsqrt(dble(ik(6)+1)/dble(ik(3)+1))
      else
        k1=1
        if(((k2/2)*2) /= k2)k1=0
        wk1=dble(k1)
        call cle0sm(bd(1),bd(3),wk1,bk(1),bk(3),a)
        if(dabs(a) < eps)return
        call rwjj(ik(3),ik(1),id(1),k1,k2,w)
        a=a*w
        wj=a/dsqrt(two*two*bk(1)+two)
        if(qm1 >= eps)return
        if(((k2/2)*2) /= k2)wj=-wj
      end if
      return
      end subroutine wj1
