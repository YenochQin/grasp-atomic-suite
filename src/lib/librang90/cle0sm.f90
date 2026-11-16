!*******************************************************************
!                                                                  *
      subroutine cle0sm(q,qm,s,c,cm,a)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   S  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                            December 1993   *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, tenth, one, two, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i
      use c1e0sm_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in)  :: q, qm, s, c, cm
      real(real64), intent(out) :: a
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iiq, iic, iis
!-----------------------------------------------
      a=zero
      iiq=two*q+tenth
      iic=two*c+tenth
      iis=two*s+tenth
      if(ittk(iiq,iic,iis).eq.0)return
      if(s.lt.eps) then
       if((q+tenth).lt.dabs(qm))return
        if((c+tenth).lt.dabs(cm))return
        if(dabs(q-c).gt.eps)return
        if(dabs(qm-cm).gt.eps)return
        a=one
      else
        call c1e0sm(q,qm,c,cm,a)
      end if
      return
      end subroutine cle0sm
