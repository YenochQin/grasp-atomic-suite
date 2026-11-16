module rmeajj_i
    implicit none (type, external)
    private
    public :: rmeajj
contains
!*******************************************************************
!                                                                  *
      subroutine rmeajj(ll,it,lq,j,its,lqs,j1s,coef)
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
      use cons_c,          only: zero
      use ribojj_c, only: imgnjj, imgtjj, impnjj, imptjj
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rmeajj11_i, only: rmeajj11
      use rmeajj9_i, only: rmeajj9
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: ll, it, lq, j, its, lqs, j1s
      real(real64), intent(out) :: coef
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i1, i2, ie1,  jt, jts
      integer :: ids3(1, 2), idv3(1, 2)
      integer :: ids5(3, 3), idv5(3, 3)
      integer :: ids7(6, 8), idv7(6, 8)
!-----------------------------------------------
      data ids3/12,20/
      data idv3/1,1/
      data ids5/-24,2*0,-30,120,-90,-54,-48,330/
      data idv5/4*1,2*7,1,2*7/
      data ids7/2*0,-40,3*0,-54,-33,-40,-65,30,0,88,-1,0,-1755,-40,  &
      0,198,-36,-72,4500,-234,-360,-78,156,0,108,-30,96,66,49,0,27,  &
      130,136,0,195,-104,-245,-624,1224,3*0,2720,-2448,-7752/
      data idv7/6*1,7,2*1,7,2*1,7,2*1,77,11,1,7,11,1,77,2*11,35,5,1, &
      91,1,13,2*5,2*1,2*7,1,11,1,3*11,3*1,143,77,91/
!
      coef=zero
      if(ll > 37) return
      if(ll <= 9) then
        if(it < 300) then
          if(imptjj(it) /= impnjj(its)) return
          if(it < its) then
            jt=it
            jts=its
          else
            jt=its
            jts=it
          end if
        end if
      end if
      if(ll == 1) then
        coef=-2
      else if(ll == 3)then
        i1=jt-2
        i2=jts-3
        coef=-dsqrt(dble(ids3(i1,i2))/dble(idv3(i1,i2)))
      else if(ll == 5)then
        i1=jt-5
        i2=jts-8
        if(ids5(i1,i2) >= 0) then
          coef=dsqrt(dble(ids5(i1,i2))/dble(idv5(i1,i2)))
        else
          coef=-dsqrt(-dble(ids5(i1,i2))/dble(idv5(i1,i2)))
        end if
      else if(ll == 7)then
        i1=jt-11
        i2=jts-17
        if(ids7(i1,i2) >= 0) then
          coef=dsqrt(dble(ids7(i1,i2))/dble(idv7(i1,i2)))
        else
          coef=-dsqrt(-dble(ids7(i1,i2))/dble(idv7(i1,i2)))
        end if
      else if(ll == 9) then
        if(it > 300) then
          if(its < 300) then
            write(6,"(A)") " ERROR IN  RMEAJJ "
            stop
          end if
          if(ll == j1s) then
            call rmeajj11(it,its,ll,coef)
          else
            call rmeajj11(its,it,ll,coef)
            ie1=lq-lqs+j-j1s+ll-1
            if((ie1/4)*4 /= ie1)coef=-coef
          end if
        else
          call rmeajj9(it,lq,j,its,lqs,j1s,coef)
          if(it > its) then
            if(mod(lq+j-lqs-j1s+ll-1,4) /= 0) coef=-coef
          end if
          write(0,"(A)") " KLAIDA RMEAJJ SUB. "
          stop
        end if
      else
        if(ll == j1s) then
          call rmeajj11(it,its,ll,coef)
        else
          call rmeajj11(its,it,ll,coef)
          ie1=lq-lqs+j-j1s+ll-1
          if((ie1/4)*4 /= ie1)coef=-coef
        end if
      end if
      if(ll < 9) then
        if(it > its) then
          if(mod(lq+j-lqs-j1s+ll-1,4) /= 0) coef=-coef
        end if
      end if
      return
      end subroutine rmeajj
end module rmeajj_i
