!***********************************************************************
!                                                                      *
      real(real64) function eigen (j)
!                                                                      *
!   This function computes an estimate of the energy of orbital J .    *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c, only: double, max_radial_points, speed_of_light_au
      use grid_c, only: double, radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orb_c, only: double, max_orbitals, orbital_kappa
      use pote_c, only: double, max_radial_points, xp, xq, yp
      use tatb_c, only: double, mtp, radial_workspace_capacity, ta, tb
      use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quad_i, only: quad
      use dpbdt_i, only: dpbdt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: piece1, piece2, piece3, piece4, piece5
!-----------------------------------------------
!
!
!   Initialization
!
      mtp = orbital_last_grid_index(j)
!
!   Exchange term
!
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = (radial_large_components(i,j)*xq(i)-radial_small_components(i,j)*xp(i))*radial_derivative_over_radius(i)
      end do
      call quad (piece1)
      piece1 = speed_of_light_au*piece1
!
!   Direct term
!
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = (radial_large_components(i,j)**2 + radial_small_components(i,j)**2)*yp(i)*radial_derivative_over_radius(i)
      end do
      call quad (piece2)
!
!   Kinetic energy terms
!
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = radial_small_components(i,j)**2*radial_grid_derivative(i)
      end do
      call quad (piece3)
      piece3 = 2.0_real64*speed_of_light_au*speed_of_light_au*piece3
!
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = (radial_large_components(i,j)*radial_small_components(i,j))*radial_derivative_over_radius(i)
      end do
      call quad (piece4)
      piece4 = -2.0_real64*dble(orbital_kappa(j))*speed_of_light_au*piece4
!
      call dpbdt (j)
      ta(1) = 0.0_real64
      do i = 2, mtp
         ta(i) = radial_large_components(i,j)*tb(i) - radial_small_components(i,j)*ta(i)
      end do
      call quad (piece5)
      piece5 = speed_of_light_au*piece5/transformed_grid_step
!
!   Assembly
!
      eigen = piece1 + piece2 + piece3 + piece4 + piece5
!
      return
!
      end function eigen
