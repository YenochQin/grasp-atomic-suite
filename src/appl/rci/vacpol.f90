!***********************************************************************
!                                                                      *
      subroutine vacpol
!                                                                      *
!   This  routine controls the setting up of the vacuum polarization   *
!   potential for the given nuclear charge distribution at each grid   *
!   point using the analytic functions defined by  L Wayne Fullerton   *
!   and G A Rinker Jr in Phys Rev A 13 (1976) 1283. The potential is   *
!   accumulated  in  array TB(I), I = 1, ..., N  which  is in COMMON   *
!   block /TATB/ .                                                     *
!                                                                      *
!   Call(s) to: [RCI92]: VAC2, VAC4.                                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use ncdist_c, only: double, max_radial_points, zdist
      use tatb_c, only: double, mtp, radial_workspace_capacity, ta, tb
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use vac2_i, only: vac2
      use vac4_i, only: vac4
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!
!   Redefine  ZDIST  to be  rho*r*r'
!
      zdist(:mtp) = zdist(:mtp)*radial_grid(:mtp)*radial_grid_derivative(:mtp)
!
!   Second-order vacuum polarisation potential; returned in
!   array TB
!
      call vac2
!
!   Fourth-order vacuum polarization potential; returned in
!   array TA
!
      call vac4
!
!   If option 7 is set, use user-defined vacuum polarization
!   potential
!
      return
      end subroutine vacpol
