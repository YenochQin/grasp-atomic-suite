module rmeajj11_i
    implicit none (type, external)
    private
    public :: rmeajj11
contains
!*******************************************************************
!                                                                  *
      subroutine rmeajj11(j1,j2,ll,s)
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
      use cons_c,          only: zero, half, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rumtjj_i, only: rumtjj
      use c0t5s_i, only: c0t5s
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: ll, j1, j2
      real(real64), intent(out) :: s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: lq, lv, lqs, lvs, j, j1s, n, inn
      real(real64) :: a1, a4, q, qq, qs, qm, qms
!-----------------------------------------------
      s=zero
      call rumtjj(j1,ll,lq,lv,j)
      call rumtjj(j2,ll,lqs,lvs,j1s)
      qq=half
      n=2
      q=half*dble(lq)
      qs=half*dble(lqs)
      qm=-half*dble(half*(ll+1)-n)
      qms=-half*dble(half*(ll+1)-(n-1))
      call c0t5s(qs,qms,qq,q,qm,a4)
      if(dabs(a4) < eps) return
      a1=dble(n*(lq+1)*(j+1))
      s=dsqrt(a1)/a4
      inn=-n-1
      if((inn/2)*2 /= inn)s=-s
      return
      end subroutine rmeajj11
end module rmeajj11_i
