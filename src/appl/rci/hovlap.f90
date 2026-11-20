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
      use grasp_parameters,   only: max_radial_points
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use hydrogenic_radial_components, only: ph, qh
      use radial_quadrature_workspace, only: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dirac_coulomb_radial_orbital, only: dcbsrw
      use radial_quadrature, only: quad
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
