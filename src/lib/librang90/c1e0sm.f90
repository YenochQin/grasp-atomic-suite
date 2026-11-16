!*******************************************************************
!                                                                  *
      subroutine c1e0sm(q,qm,c,cm,a)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use cons_c,          only: zero, tenth, one, two, three, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), intent(in)  :: q, qm, c, cm
      real(double), intent(out) :: a
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iiq, iic, is, ig
!-----------------------------------------------
      a=zero
      iiq=two*q+tenth
      iic=two*c+tenth
      if(ittk(iiq,iic,2) == 0)return
      if(dabs(qm-cm) > eps) return
      if((q+tenth) < dabs(qm)) return
      if((c+tenth) < dabs(cm)) return
      if(dabs(qm) <= eps) then
       is=q+c+one+tenth
       if((is/2)*2 /= is) return
      end if
      ig=q-c+two+tenth
      if(ig <= 0) return
      if(ig > 3) return
      if (ig == 1) then
        a=dsqrt(((c+cm)*(c-cm))/((two*c-one)*c))
      else if (ig == 2) then
        a=cm/dsqrt(c*(c+one))
      else if (ig == 3) then
        a=-dsqrt(((c+cm+one)*(c-cm+one))/((c+one)*(two*c+three)))
      end if
      return
      end subroutine c1e0sm
