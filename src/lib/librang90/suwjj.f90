!*******************************************************************
!                                                                  *
      subroutine suwjj(k1,k2,ll,j1,j2,suw)
!                                                                  *
!                 (k1 k2)                                          *
!     ( j QJ ::: W      ::: j QJ )                                 *
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
      use cons_c,          only: zero, eps
      use ribojj_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rumtjj_i
      use ixjtik_i
      use rmeajj_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: k1, k2, j1, j2
      integer,      intent(out) :: ll
      real(real64), intent(out) :: suw
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kk1,kk2,i,ip,ig,l2q1,l2v1,l2j1,l2q2,l2v2,l2j2, &
                 l2qi,l2vi,l2ji
      real(real64) :: coef1, coef2, s, si1, si2
!-----------------------------------------------
      suw=zero
      if(imptjj(j1) /= imptjj(j2)) return
      s=zero
      call rumtjj(j1,ll,l2q1,l2v1,l2j1)
      kk1=k1*2
      kk2=k2*2
      ip=impnjj(j1)
      ig=imgnjj(j1)
      call rumtjj(j2,ll,l2q2,l2v2,l2j2)
      do i=ip,ig
        call rumtjj(i,ll,l2qi,l2vi,l2ji)
        if(ixjtik(ll,ll,kk2,l2j2,l2j1,l2ji) /= 0) then
          if(ixjtik(1,1,kk1,l2q2,l2q1,l2qi) /= 0) then
            call rmeajj(ll,j1,l2q1,l2j1,i,l2qi,l2ji,coef1)
            if(dabs(coef1) > eps) then
              call rmeajj(ll,i,l2qi,l2ji,j2,l2q2,l2j2,coef2)
              if(dabs(coef2) > eps) then
                call sixj(ll,ll,kk2,l2j2,l2j1,l2ji,0,si1)
                call sixj(1,1,kk1,l2q2,l2q1,l2qi,0,si2)
                s=s+si1*si2*coef1*coef2
              endif
            endif
          endif
        endif
      end do
      suw=s*dsqrt(dble((kk1+1)*(kk2+1)))
      if(mod(l2q1+l2j1+l2q2+l2j2+kk1+kk2,4) /= 0)suw=-suw
      return
      end subroutine suwjj
