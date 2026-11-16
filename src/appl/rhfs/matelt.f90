!***********************************************************************
!                                                                      *
      subroutine matelt(i1, k, i2, apart, gjpart, dgjpart)
!                                                                      *
!   This  routine computes  the angular part  of the  reduced matrix   *
!   elements of the magnetic and electric multipole operators as       *
!   well as for th operators associated with the g_j factor            *
!                                                                      *
!   Calls to: [LIB92]: CLRX                                            *
!                                                                      *
!   Written by Per Jonsson                Last revision: 22 Oct 1999   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:03   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/01/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: max_orbitals
      use orb_c,           only: orbital_principal_n, orbital_kappa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use clrx_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i1
      integer  :: k
      integer, intent(in) :: i2
      real(double), intent(out) :: apart
      real(double), intent(out) :: gjpart
      real(double), intent(out) :: dgjpart
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: kap1, kap2, l1, l2
      real(double) :: fase, ovlfac
!-----------------------------------------------
!
!
!   Set KAP1 and KAP2
!
      kap1 = orbital_kappa(i1)
!
      if (mod(k,2) == 1) then
         kap2 = -orbital_kappa(i2)
      else
         kap2 = orbital_kappa(i2)
      endif
!
!   Determine the l quantum numbers
!
      if (kap1 > 0) then
         l1 = kap1
      else
         l1 = (-kap1) - 1
      endif
!
      if (kap2 > 0) then
         l2 = kap2
      else
         l2 = (-kap2) - 1
      endif
!
      if (mod(l1 + k + l2,2) == 0) then
!
!   Parity selection rule satisfied
!
!   Determine the phase factor
!
         if (mod(kap1 + 1,2) == 0) then
            fase = 1.0d00
         else
            fase = -1.0d00
         endif
!
!   The other factor is \sqrt (2 j_2 + 1); since j = | \kappa | - 1/2,
!   we have 2 j_2 + 1 = 2 | \kappa |; the factor \sqrt (2 j_2 + 1)
!   has been accounted for in MCT
!
         ovlfac = fase*sqrt(dble(2*abs(kap2)))
!
         if (mod(k,2) == 1) then
!
!   These are for the magnetic multipole moments and the two operators of the g_j factor
!
            apart = (kap1 + orbital_kappa(i2))*clrx(kap1,k,kap2)*ovlfac
!
            gjpart = apart
            dgjpart = -(kap1 + orbital_kappa(i2)-1.0d00)*clrx(kap1,k,kap2)*ovlfac
         else
!
!   These are for the electric multipole moments
!
            apart = clrx(kap1,k,kap2)*ovlfac
!
         endif
!
      else
!
         apart = 0.0d00
         gjpart = 0.0d00
         dgjpart = 0.0d00
!
      endif
!
      return
      end subroutine matelt
