!*******************************************************************
!                                                                  *
      subroutine oneparticlejj1(ns,ka,jja,jjb,ja,jb,coeff)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
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
      use cons_c,          only: zero, half, eps
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use recop00_i
      use recop1_i
      use recop2_i
      use wj1_i
      use perko2_i, only: perko2
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ns,ka,jja,jjb,ja,jb
      real(double), intent(out) :: coeff
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: iat
      real(double) :: rec,wj,qm1,qm2
!-----------------------------------------------
!
!     THE CASE 11   + -
!
      coeff=zero
      if(ja == jb) then
        if(jja /= jjb) then
          call recop00(ns,ja,ja,ka,iat)
          if(iat == 0)return
        end if
        call recop1(ns,ja,ka,0,iat,rec)
        if(iat == 0)return
        call perko2(ja,ja,ja,ja,1)
        qm1=half
        qm2=-half
        call wj1(ik1,bk1,id1,bd1,ka,qm1,qm2,wj)
        if(dabs(wj) > eps) then
           call recop1(ns,ja,ka,1,iat,rec)
           coeff=wj*rec*dsqrt(dble(id1(3)+1))/dsqrt(dble(2*ka+1))
           coeff=-coeff
        end if
      end if
      return
      end subroutine oneparticlejj1
