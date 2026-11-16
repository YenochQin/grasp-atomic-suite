module c0t5s_i
    implicit none (type, external)
    private
    public :: c0t5s
contains
!*******************************************************************
!                                                                  *
      subroutine c0t5s(q,qm,sm,c,cm,a)
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                ---          ---  *
!                                                I  Q  1/2  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:              I              I  *
!                                                I  QM  SM  CM  I  *
!                                                ---          ---  *
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
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, tenth , half, one, two, eps
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i, only: ittk
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(real64), intent(in)  :: q, qm, sm, c, cm
      real(real64), intent(out) :: a
!      DIMENSION GC(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                    :: iiq, iic, ie
      real(real64) :: gc(2)
!-----------------------------------------------
      gc(1)=one
      gc(2)=-one
      a=zero
      iiq=two*q+tenth
      iic=two*c+tenth
      if(ittk(iiq,iic,1) == 0)return
      if(dabs(qm+sm-cm) > eps)return
      if((half+tenth) < dabs(sm))return
      if((q+tenth) < dabs(qm))return
      if((c+tenth) < dabs(cm))return
      ie=dabs(half-sm)+one+tenth
      if(dabs(q+half-c) < eps) then
        a=dsqrt((c+gc(ie)*cm)/(two*c))
      else
        if(dabs(q-half-c) > eps)return
        a=-gc(ie)*dsqrt((c-gc(ie)*cm+one)/(two*c+two))
      end if
      return
      end subroutine c0t5s
end module c0t5s_i
