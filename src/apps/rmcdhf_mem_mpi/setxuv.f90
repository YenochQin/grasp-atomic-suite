!***********************************************************************
!                                                                      *
      subroutine setxuv(orbital_index)
!                                                                      *
!   This  SUBROUTINE  sets  up the arrays XU and XV, for use by  the   *
!   subprograms IN and  OUT.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: dp, dq, max_radial_points
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
      use scf_radial_potentials, only: max_radial_points, xp, xq, yp
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: orbital_index
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, i
      real(real64) :: dmhh
!-----------------------------------------------
!
!   Define constants
!
      dmhh = -transformed_grid_step*0.5_real64
!
!   Set up arrays XU and XV; since XU(1), XV(1) are never used,
!   set them to some arbitrary value
!
      nm1 = radial_point_count - 1
      xu(1) = 0.0_real64
      xv(1) = 0.0_real64
      xu(2:nm1) = dmhh*(xp(3:nm1+1)*radial_derivative_over_radius(3:nm1+1)+xp(2:nm1)*radial_derivative_over_radius(2:nm1)) + dp(2&
         :nm1)
      xv(2:nm1) = dmhh*(xq(3:nm1+1)*radial_derivative_over_radius(3:nm1+1)+xq(2:nm1)*radial_derivative_over_radius(2:nm1)) + dq(2&
         :nm1)
!
      xu(radial_point_count) = dmhh*xp(radial_point_count)*radial_derivative_over_radius(radial_point_count) + &
         & dp(radial_point_count)
      xv(radial_point_count) = dmhh*xp(radial_point_count)*radial_derivative_over_radius(radial_point_count) + &
         & dq(radial_point_count)
!
      return
      end subroutine setxuv
