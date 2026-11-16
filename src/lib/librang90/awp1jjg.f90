module awp1jjg_i
    implicit none (type, external)
    private
    public :: awp1jjg
contains
!*******************************************************************
!                                                                  *
      subroutine awp1jjg(k1,bk2,qm1,qm2,qm3,ik,bk,id,bd,aw)
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
      use cons_c,          only: zero, two, tenth, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use itjj2_i
      use ixjtik_i
      use izas1_i
      use rumtjj_i
      use a1jj_i
      use w1jjg_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: bk2, qm1, qm2, qm3
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: aw
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk1, kk2, ie, iq, iq2, iq3, iqm, it, itp, itg, ibtt
      integer :: ibt(7)
      real(real64)               :: enqp, d1, si1, w
      real(real64) :: bt(3)
!-----------------------------------------------
      aw=zero
      if(izas1(id(7),bd(3),ik(7),bk(3)) == 0)return
      enqp=zero
      kk1=k1*2
      kk2=bk2+bk2+tenth*bk2
      iq2=qm2*two+qm2*tenth
      iq3=qm3*two+qm3*tenth
      iq=iq2+iq2
      if(id(3) > 37) return
      if(itjj2(ik,id,kk2,bk,bd,ibt,bt,itp,itg,iq) == 0)return
      iqm=two*dabs(bt(3))+tenth
      do it=itp,itg
        call rumtjj(it,ibt(3),ibt(7),ibtt,ibt(6))
        if(iqm <= ibt(7)) then
          if(ixjtik(ik(3),kk1,kk2,id(6),ik(6),ibt(6)) /= 0) then
            ibt(1)=it
            bt(2)=dble(ibt(6))/two
            bt(1)=dble(ibt(7))/two
            call a1jj(ik,bk,ibt,bt,qm1,d1)
            if(dabs(d1) > eps) then
              call w1jjg(k1,qm2,qm3,ibt,bt,id,bd,w)
              if(dabs(w) > eps) then
                d1=d1*w
                call sixj(ik(3),kk1,kk2,id(6),ik(6),ibt(6),0,si1)
                enqp=enqp+d1*si1
              endif
            endif
          endif
        endif
      end do
      aw=enqp*dsqrt(dble((kk2+1)))
      ie=kk2+ik(6)+id(6)
      if(((ie/4)*4) /= ie)aw=-aw
      return
      end subroutine awp1jjg
end module awp1jjg_i
