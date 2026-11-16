module w1jjg_i
    implicit none (type, external)
    private
    public :: w1jjg
contains
!*******************************************************************
!                                                                  *
      subroutine w1jjg(k1,qm1,qm2,ik,bk,id,bd,ww)
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
      use cons_c,          only: zero, tenth, half, two, eps
      use ribojj_c
      use ribojj9_c
      use ribojj11_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use rumtjj_i
      use ixjtik_i
      use ittk_i
      use a1jj_i
      use sixj_i
      use izas1_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)               :: k1
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in)               :: qm1, qm2
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out)              :: ww
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk1,iq,iqm,ie,it,itk,itd,itp,itp1,itg,itg1,ibtt
      integer :: ibt(7)
      real(real64)               :: enqp, d1, si1, w
      real(real64) :: bt(3)
!-----------------------------------------------
      ww=zero
      if(izas1(id(7),bd(3),ik(7),bk(3)) == 0)return
      enqp=zero
      kk1=k1*2
      iq=qm2*two+qm2*tenth
      if(id(3) > 37) return
      if(ittk(id(6),ik(6),kk1) == 0)return
      itk=ik(1)
      itd=id(1)
      if(id(3) == 9) then
        if(id(4) > 2) call mes(1)
        if(ik(4) > 2) call mes(1)
        itk=itk-300
        itd=itd-300
        itp1=impnjj9(itk)
        itp=impnjj9(itd)
        if(itp1 /= itp)return
        itg1=imgnjj9(itk)
        itg=imgnjj9(itd)
       else
        if(id(4) > 2) call mes(1)
        if(ik(4) > 2) call mes(1)
        itp1=impnjj11(itk)
        itp=impnjj11(itd)
        if(itp1 /= itp)return
        itg1=imgnjj11(itk)
        itg=imgnjj11(itd)
      endif
      if(itg1 /= itg)return
      ibt(2)=id(2)
      ibt(3)=id(3)
      ibt(4)=id(4)+iq
      bt(3)=bd(3)+half*dble(iq)
      iqm=two*dabs(bt(3))+tenth
      do it=itp,itg
        call rumtjj(it,ibt(3),ibt(7),ibtt,ibt(6))
        if(iqm <= ibt(7)) then
          if(ixjtik(ik(3),ik(3),kk1,id(6),ik(6),ibt(6)) /= 0) then
            ibt(1)=it
            bt(2)=dble(ibt(6))/two
            bt(1)=dble(ibt(7))/two
            call a1jj(ik,bk,ibt,bt,qm1,d1)
            if(dabs(d1) > eps) then
              call a1jj(ibt,bt,id,bd,qm2,w)
              if(dabs(w) > eps) then
                d1=d1*w
                call sixj(ik(3),ik(3),kk1,id(6),ik(6),ibt(6),0,si1)
                enqp=enqp+d1*si1
              endif
            endif
          endif
        endif
      end do
      ww=enqp*dsqrt(dble((kk1+1)))
      ie=kk1+ik(6)+id(6)
      if(((ie/4)*4) /= ie)ww=-ww
      return
      end subroutine w1jjg
end module w1jjg_i
