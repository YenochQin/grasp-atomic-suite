module a1jj_i
    implicit none (type, external)
    private
    public :: a1jj
contains
!*******************************************************************
!                                                                  *
      subroutine a1jj(ik,bk,id,bd,qm1,a)
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
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in) :: ik(7), id(7)
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(in)               :: qm1
      real(real64), intent(out)              :: a
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ifaz, isuma
      real(real64) :: ab
!-----------------------------------------------
      a=zero
      if(qm1 < eps) then
        isuma=(id(6)+1)*id(4)
        ab=dble(isuma)
        a=dsqrt(ab)
        ifaz=id(6)+id(3)-ik(6)+id(4)*2
        if((ifaz/4)*4 /= ifaz)a=-a
      else
        isuma=(ik(6)+1)*ik(4)
        ab=dble(isuma)
        a=dsqrt(ab)
        if((ik(4)/2)*2 /= ik(4))a=-a
      end if
      return
      end subroutine a1jj
end module a1jj_i
