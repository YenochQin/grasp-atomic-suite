!***********************************************************************
!                                                                      *
      subroutine setxv(j)
!                                                                      *
!   This  subprogram  sets up the inhomogeneous terms for the varia-   *
!   tion equations.                                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: max_radial_points, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: hhc
!-----------------------------------------------
!
      hhc = 0.5_real64*transformed_grid_step/speed_of_light_au
!
!   Set up arrays TF and TG
!
      do i = 1, radial_point_count
         xu(i) = -radial_small_components(i,j)*hhc*radial_grid_derivative(i)
         xv(i) = radial_large_components(i,j)*hhc*radial_grid_derivative(i)
      end do
!
      xu(:radial_point_count-1) = xu(2:radial_point_count) + xu(:radial_point_count-1)
      xv(:radial_point_count-1) = xv(2:radial_point_count) + xv(:radial_point_count-1)
!
      return
!
      end subroutine setxv
