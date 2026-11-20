!***********************************************************************
!                                                                      *
      subroutine dampck(ipr, orbital_index, ed1, ed2)
!                                                                      *
!   This subroutine determines the damping factor appropriate to the   *
!   present  orbital. The algorithm is taken from C Froese Fischer's   *
!   program MCHF.                                                      *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
      use orbital_csf_basis, only: max_orbitals, orbital_energy_parameters, ped
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
      ed2 = (ed2-orbital_energy_parameters(orbital_index))/ed2
      if (adaptv) then
         if (ed1*ed2 < -1.0e-4_real64) then
!           We have a significant oscillation, damp
            orbital_damping_factors(orbital_index) = 0.10_real64 + 0.90_real64*orbital_damping_factors(orbital_index)
!           Print *, 'Increasing ODAMP', odamp(j)
         else
            orbital_damping_factors(orbital_index) = 0.50_real64*orbital_damping_factors(orbital_index)
         end if
      end if
!
!     Save the relative difference for the next update of this orbital
      ped(orbital_index) = ed2
!
      return
      end subroutine dampck
