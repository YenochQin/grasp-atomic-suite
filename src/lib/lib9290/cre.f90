!***********************************************************************
!                                                                      *
      real(real64) function cre (kap1, k, kap2)
!-----------------------------------------------
!                                                                      *
!   Computes the relativistic reduced matrix element                   *
!                                                                      *
!                         (j1 || C(K) || j2),                          *
!                                                                      *
!   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
!   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
!   conditions are tested by the routine CLRX.                         *
!                                                                      *
!   Call(s) to: [LIB92] CLRX.                                          *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use clrx_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: kap1
      integer, intent(in) :: k
      integer, intent(in) :: kap2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k1
      real(real64) :: dk1k2
!-----------------------------------------------
!
      k1 = abs(kap1)
      dk1k2 = dble(4*k1*iabs(kap2))
      cre = sqrt(dk1k2)*clrx(kap1,k,kap2)
      if (mod(k1,2) == 1) cre = -cre
!
      return
      end function cre
