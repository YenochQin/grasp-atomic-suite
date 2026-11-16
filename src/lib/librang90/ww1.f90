!*******************************************************************
!                                                                  *
      subroutine ww1(ik,bk,id,bd,k2,qm1,qm2,qm3,qm4,ww)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 25  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                    N      (k2)   (k2) (0)   N'     +-            *
!     ELEMENT      (j QJ::[W   *  W    ]   ::j QJ)   -+            *
!                                                    ++            *
!                                                    -- B17 (2.4)  *
!                                                                  *
!     SUBROUTINE CALLED: ITJJ,IXJTIK,IZAS1,RUMTJJ,WJ1              *
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
      use cons_c,          only: zero, tenth, two, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use izas1_i
      use itjj_i
      use rumtjj_i
      use ixjtik_i
      use wj1_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k2
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in)               :: qm1, qm2, qm3, qm4
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(out)              :: ww
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk2,kk6,iq,iqm,iq3,iq4,ie,ie1,it,itp,itg,ibtt
      integer, dimension(7)      :: ibt
      real(double)               :: enqp, d1, w
      real(double), dimension(3) :: bt
!-----------------------------------------------
      ww=zero
      if(id(6) /= ik(6))return
      if(izas1(id(7),bd(3),ik(7),bk(3)) == 0)return
      enqp=zero
      kk2=k2*2
      iq3=qm3*two+qm3*tenth
      iq4=qm4*two+qm4*tenth
      iq=iq3+iq4
      if(itjj(ik,id,0,bk,bd,ibt,bt,kk6,itp,itg,iq) == 0)return
      ie1=kk2-ik(6)
      iqm=two*dabs(bt(3))+tenth
      do it=itp,itg
        call rumtjj(it,ibt(3),ibt(7),ibtt,ibt(6))
        if(iqm > ibt(7)) cycle
        if(ixjtik(kk2,kk2,0,id(6),ik(6),ibt(6)) == 0) cycle
        ibt(1)=it
        bt(2)=dble(ibt(6))/two
        bt(1)=dble(ibt(7))/two
        call wj1(ik,bk,ibt,bt,k2,qm1,qm2,d1)
        if(dabs(d1) < eps) cycle
        call wj1(ibt,bt,id,bd,k2,qm3,qm4,w)
        if(dabs(w) < eps) cycle
        d1=d1*w
        ie=ie1+ibt(6)
        if(((ie/4)*4) /= ie)d1=-d1
        enqp=enqp+d1
      end do
      ww=enqp/dsqrt(dble(kk2+1)*(ik(6)+1))
      return
      end subroutine ww1
