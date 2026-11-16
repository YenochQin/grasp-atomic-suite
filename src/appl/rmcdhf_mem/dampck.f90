!***********************************************************************
!                                                                      *
      subroutine dampck(ipr, j, ed1, ed2)
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
      use damp_c
      use orb_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(inout) :: ipr
      integer, intent(in) :: j
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
      adaptv = orbital_damping_factors(j) >= 0.0_real64
!
      ed2 = (ed2-orbital_energy_parameters(j))/ed2
      if (adaptv) then
         if (ed1*ed2 < -1.0e-4_real64) then
!           We have a significant oscillation, damp
            orbital_damping_factors(j) = 0.10_real64 + 0.90_real64*orbital_damping_factors(j)
!           Print *, 'Increasing ODAMP', odamp(j)
         else
            orbital_damping_factors(j) = 0.50_real64*orbital_damping_factors(j)
         endif
      endif
!
!     Save the relative difference for the next update of this orbital
      ped(j) = ed2
!
      return
      end subroutine dampck
