!*******************************************************************
!                                                                  *
      subroutine c1e1sm(q,qm,sm,c,cm,a)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  1  CM  I  *
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
      real(double), intent(in)  :: q, qm, sm, c, cm
      real(double), intent(out) :: a
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                    :: ie, iiq, iic
      real(double), dimension(2) :: gc
!-----------------------------------------------
      gc(1)=one
      gc(2)=-one
      a=zero
      iiq=two*q+tenth
      iic=two*c+tenth
      if(ittk(iiq,iic,2).eq.0)return
      if(dabs(qm+sm-cm).gt.eps)return
      if((q+tenth).lt.dabs(qm))return
      if((c+tenth).lt.dabs(cm))return
      ie=0
      if(dabs(sm-one).lt.eps)ie=1
      if(dabs(sm+one).lt.eps)ie=2
      if(ie.eq.0)return
      if(dabs(q+one-c).lt.eps) then
        a=dsqrt((c+gc(ie)*cm-one)*(c+gc(ie)*cm)/ &
        ((two*c-one)*two*c))
      else if(dabs(q-c).lt.eps) then
        a=-gc(ie)*dsqrt((c+gc(ie)*cm)*(c-gc(ie)*cm+one)/ &
        ((c+one)*two*c))
      else if(dabs(q-one-c).gt.eps) then
        return
      else
        a=dsqrt((c-gc(ie)*cm+one)*(c-gc(ie)*cm+two)/ &
        ((two*c+two)*(two*c+three)))
      end if
      return
      end subroutine c1e1sm
