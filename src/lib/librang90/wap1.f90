!*******************************************************************
!                                                                  *
      subroutine wap1(ik,bk,id,bd,k1,bk2,qm1,qm2,qm3,wa)
!                                                                  *
!   ---------------  SECTION SQ    SUBPROGRAM 23  --------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
!                                                                  *
!                    N       (k1)   (j) (k2)   N'    +-            *
!     ELEMENT:     (j QJ::[ W    * A   ]    ::j QJ)  -+            *
!                                                    ++            *
!                                                    -- B17 (2.3)  *
!                                                                  *
!     SUBROUTINE CALLED: C0T5S,ITJJ3,IXJTIK,IZAS1,RUMTJJ,SIXJ,     *
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
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, tenth, two, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use wap1jjg_i
      use izas1_i
      use itjj3_i
      use rumtjj_i
      use ixjtik_i
      use c0t5s_i
      use rmeajj_i
      use wj1_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: bk2, qm1, qm2, qm3
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk1, kk2, ie, iq3, iqm, it, itp, itg, ibtt
      integer :: ibt(7)
      real(real64)               :: enqp, d1, s, si, w
      real(real64) :: bt(3)
!-----------------------------------------------
      wa=zero
      if(id(3) == 9) then
        if(max0(ik(4),id(4)) < 3) then
          if(ik(1) < 300) call mes(55)
          if(id(1) < 300) call mes(55)
          call wap1jjg(k1,bk2,qm1,qm2,qm3,ik,bk,id,bd,wa)
          return
        endif
      elseif(id(3) > 9) then
        call wap1jjg(k1,bk2,qm1,qm2,qm3,ik,bk,id,bd,wa)
        return
      endif
      if(izas1(id(7),bd(3),ik(7),bk(3)) == 0)return
      enqp=zero
      kk1=k1*2
      iq3=qm3*two
      kk2=bk2+bk2+tenth*bk2
      if(itjj3(ik,id,kk2,bk,bd,ibt,bt,itp,itg,iq3) == 0)return
      iqm=two*dabs(bt(3))+tenth
      do it=itp,itg
        call rumtjj(it,ibt(3),ibt(7),ibtt,ibt(6))
        if(iqm > ibt(7)) cycle
        if(ixjtik(kk1,ik(3),kk2,id(6),ik(6),ibt(6)) == 0) cycle
        ibt(1)=it
        bt(2)=dble(ibt(6))/two
        bt(1)=dble(ibt(7))/two
        call c0t5s(bd(1),bd(3),qm3,bt(1),bt(3),d1)
        if(dabs(d1) < eps) cycle
        call rmeajj(id(3),ibt(1),ibt(7),ibt(6),id(1),id(7),id(6),s)
        if(dabs(s) < eps) cycle
        d1=d1*s
        call wj1(ik,bk,ibt,bt,k1,qm1,qm2,w)
        if(dabs(w) < eps) cycle
        d1=d1*w
        call sixj(kk1,ik(3),kk2,id(6),ik(6),ibt(6),0,si)
        d1=d1*si/dsqrt(dble(ibt(7)+1))
        enqp=enqp+d1
      end do
      wa=enqp*dsqrt(dble(kk2+1))
      ie=kk2+ik(6)+id(6)+2
      if(((ie/4)*4) /= ie)wa=-wa
      return
      end subroutine wap1
