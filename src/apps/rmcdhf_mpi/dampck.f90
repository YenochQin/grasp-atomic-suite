!***********************************************************************
!                                                                      *
      subroutine dampck(ipr, orbital_index, ed1, ed2)
!                                                                      *
!   This subroutine determines the damping factor appropriate to the   *
!   present  orbital. The algorithm is taken from C Froese Fischer's   *
!   program MCHF.                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!   Modified by C. Froese Fischer           Last update: 07 Apr 2009   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
      use orbital_csf_basis, only: max_orbitals, orbital_energy_parameters
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: ipr
      integer, intent(in) :: orbital_index
      real(real64), intent(inout) :: ed1
      real(real64), intent(inout) :: ed2
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: adaptv
!-----------------------------------------------
!
!   The damping is adaptive (i.e., can be modified by this SUBROUTINE)
!   if and only if ODAMP(J) .GE. 0.0
!
      adaptv = orbital_damping_factors(orbital_index) >= 0.0_real64
!
      if (ipr /= orbital_index) then
         if (adaptv) orbital_damping_factors(orbital_index) = 0.75_real64*orbital_damping_factors(orbital_index)
      else
         ed2 = ed2-orbital_energy_parameters(orbital_index)
         if (adaptv) then
            if (ed1*ed2 > 0.0_real64) then
               orbital_damping_factors(orbital_index) = 0.75_real64*orbital_damping_factors(orbital_index)
            else
               orbital_damping_factors(orbital_index) = 0.25_real64+0.75_real64*orbital_damping_factors(orbital_index)
            end if
         end if
      end if
!
      if (ipr == orbital_index) then
         ed1 = ed2
      else
         ed1 = ed2-orbital_energy_parameters(orbital_index)
      end if
!
      ipr = orbital_index
      return
      end subroutine dampck
