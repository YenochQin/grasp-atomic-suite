!*******************************************************************
!                                                                  *
      subroutine awp1(ik,bk,id,bd,k1,bk2,qm1,qm2,qm3,aw)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 01  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                    N      (j)  (k1) (k2)   N'     +-             *
!     ELEMENT:     (j QJ ::[A  * W    ]   ::j  QJ)  -+             *
!                                                   ++             *
!                                                   --  B17 (2.3)  *
!                                                                  *
!     SUBROUTINE CALLED: C0T5S,ITJJ2,IXJTIK,IZAS1,RUMTJJ,SIXJ,     *
!                        RMEAJJ,WJ1                                *
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
      use cons_c,          only: zero, two, tenth, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use itjj2_i
      use ixjtik_i
      use izas1_i
      use awp1jjg_i, only: awp1jjg
      use c0t5s_i, only: c0t5s
      use mes_i, only: mes
      use rmeajj_i, only: rmeajj
      use rumtjj_i, only: rumtjj
      use sixj_i, only: sixj
      use wj1_i, only: wj1
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in)               :: bk2, qm1, qm2, qm3
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(out)              :: aw
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk1, kk2, ie, iq, iq2, iq3, iqm, it, itp, itg, ibtt
      integer,      dimension(7) :: ibt
      real(double)               :: enqp, d1, s, si, w
      real(double), dimension(3) :: bt
!-----------------------------------------------
      aw=zero
      if(id(3) == 9) then
        if(max0(ik(4),id(4)) < 3) then
          if(ik(1) < 300) call mes(54)
          if(id(1) < 300) call mes(54)
          call awp1jjg(k1,bk2,qm1,qm2,qm3,ik,bk,id,bd,aw)
          return
        else
          print*, "ERROR in AWP1"
          stop
        endif
      elseif(id(3) > 9) then
        call awp1jjg(k1,bk2,qm1,qm2,qm3,ik,bk,id,bd,aw)
        return
      endif
      if(izas1(id(7),bd(3),ik(7),bk(3)) == 0)return
      enqp=zero
      iq2=qm2*two+qm2*tenth
      iq3=qm3*two+qm3*tenth
      iq=iq2+iq3
      kk1=k1*2
      kk2=bk2+bk2+tenth*bk2
      if(itjj2(ik,id,kk2,bk,bd,ibt,bt,itp,itg,iq) == 0)return
      iqm=two*dabs(bt(3))+tenth
      do it=itp,itg
        call rumtjj(it,ibt(3),ibt(7),ibtt,ibt(6))
        if(iqm > ibt(7)) cycle
        if(ixjtik(ik(3),kk1,kk2,id(6),ik(6),ibt(6)) == 0) cycle
        ibt(1)=it
        bt(2)=dble(ibt(6))/two
        bt(1)=dble(ibt(7))/two
        call c0t5s(bt(1),bt(3),qm1,bk(1),bk(3),d1)
        if(dabs(d1) < eps) cycle
        call rmeajj(ik(3),ik(1),ik(7),ik(6),ibt(1),ibt(7),ibt(6),s)
        if(dabs(s) < eps) cycle
        call wj1(ibt,bt,id,bd,k1,qm2,qm3,w)
        d1=d1*w*s
        if(dabs(d1) < eps) cycle
        call sixj(ik(3),kk1,kk2,id(6),ik(6),ibt(6),0,si)
        d1=d1*si/dsqrt(dble(ik(7)+1))
        enqp=enqp+d1
      end do
      aw=enqp*dsqrt(dble(kk2+1))
      ie=kk2+ik(6)+id(6)+2
      if(((ie/4)*4) /= ie)aw=-aw
      return
      end subroutine awp1
