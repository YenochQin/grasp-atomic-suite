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
      use vast_kind_param, only: double
      use cons_c,          only: zero, eps
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in), dimension(7) :: ik, id
      real(double), intent(in), dimension(3) :: bk, bd
      real(double), intent(in)               :: qm1
      real(double), intent(out)              :: a
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ifaz, isuma
      real(double) :: ab
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
      endif
      return
      end subroutine a1jj
