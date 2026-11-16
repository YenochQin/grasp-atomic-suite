!***********************************************************************
!                                                                      *
      real(real64) function hovlap (p, q, mtpo, np, kappa, z)
!                                                                      *
!   This subprogram computes the overlap of the orbital tabulated in   *
!   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
!   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
!                                                                      *
!   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use grid_c
      use horb_c, only: ph, qh
      use tatb_c, only: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dcbsrw_i, only: dcbsrw
      use quad_i, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mtpo
      integer  :: np
      integer  :: kappa
      real(real64)  :: z
      real(real64), intent(in) :: p(max_radial_points)
      real(real64), intent(in) :: q(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtph, i
      real(real64) :: eh, pzh, result
!-----------------------------------------------
!
!   Set up the hydrogenic orbital
!
      call dcbsrw (np, kappa, z, eh, pzh, ph, qh, mtph)
!
!   Compute the overlap
!
      mtp = min(mtph,mtpo)
      ta(1) = 0.0_real64
      ta(2:mtp) = (p(2:mtp)*ph(2:mtp)+q(2:mtp)*qh(2:mtp))*radial_grid_derivative(2:mtp)
      call quad (result)
!
      hovlap = result
!
      return
      end function hovlap
