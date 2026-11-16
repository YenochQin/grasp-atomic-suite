!***********************************************************************
!                                                                      *
      real(real64) function shield (j)
      use, intrinsic :: iso_fortran_env, only: real32
!                                                                      *
!   This  routine estimates the screening (or shielding) for orbital   *
!   J according to the relativistic hydrogenic energy formula.         *
!                                                                      *
!   Written by Yu Zou, at Vanderbilt University: 02 Mar 2000           *
!                                                                      *
!   Formula:                                                           *
!    epsilon = E(n,k)/Z^2_{eff}                                        *
!    beta = alpha Z_{eff}                                              *
!    epsilon = 1/beta^{2} [(1+beta^2/nu^2)^{-1/2} - 1]                 *
!    nu = n + (k^2 - beta^2)^{1/2} - |k|                               *
!   where E(n,k) is the orbital energy for the specific principal      *
!   quantum n and k quantum number, alpha is the fine structure        *
!   constant.                                                          *
!   Note that SHIELD is set to -C if Zeff is out of range. This will   *
!   ensure that the self-energy for such orbital is zero.              *
!   Reference: M J Seaton, Rep. Prog. Phys., Vol.46, P167,1983         *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c
      use orb_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: a, aa, bb, cc, t, beta
!-----------------------------------------------
!
!
! a=[(1+\alpha^2 E)^{-2} - 1]^{-1/2}
! note E(J) = -E(n,k)
      a = 1 - orbital_energy_parameters(j)/(speed_of_light_au*speed_of_light_au)
      a = 1/(a*a) - 1
      if (a < 0.0_real32) then
         shield = -speed_of_light_au
         return
      endif
      a = 1/sqrt(a)
!
      aa = 1 + a*a
      bb = a*dble(orbital_principal_n(j)-iabs(orbital_kappa(j)))
      cc = dble((orbital_principal_n(j)-iabs(orbital_kappa(j)))**2-orbital_kappa(j)**2)
      t = bb*bb - aa*cc
      if (t < 0.0_real32) then
         shield = -speed_of_light_au
         return
      endif
      beta = bb + sqrt(t)
      beta = beta/aa
      shield = nuclear_charge - beta*speed_of_light_au
!     print *, j, '  Initial Screen Factor = ',shield
      return
      end function shield
