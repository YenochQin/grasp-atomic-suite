!***********************************************************************
!                                                                      *
      subroutine prwf(j)
!                                                                      *
!   Makes a (debug) printout of wave functions. There are two modes:   *
!                                                                      *
!      (1) J>0  - Used as a debug option in SOLVE, wavefunctions for   *
!                 orbital J are printed                                *
!      (2) J=0  - A  printout  of the grid and all wave functions is   *
!                 made                                                 *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      use atomic_calculation_state, only: max_radial_points, nuclear_charge, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
      use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_principal_n, orbital_symmetry_labels
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE draw_I
      use radial_ascii_plot, only: draw
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nb2, nrows, ii, ii1, ii2, k, mfk, jgg
      real(real64) :: cbz
!-----------------------------------------------
!
!
      cbz = speed_of_light_au/nuclear_charge
!
      if (j > 0) then
!
!   Mode (1)
!
         write (99, '("1",1i2,1a2," orbital:"/,2(" --------- r --------- ------- P (r) -------",&
            &" ------- Q (r) -------"))') orbital_principal_n(j), orbital_symmetry_labels(j)
         nb2 = mtp0/2
         if (2*nb2 == mtp0) then
            nrows = nb2
         else
            nrows = nb2 + 1
         end if
         do ii = 1, nrows
            ii1 = ii
            ii2 = ii1 + nrows
            if (ii2 <= mtp0) then
               write (99, '(1p,6(1x,1d21.14))') radial_grid(ii1), p(ii1), q(ii1), radial_grid(ii2), p(ii2), q(ii2)
            else if (ii1 <= mtp0) then
               write (99, '(1p,6(1x,1d21.14))') radial_grid(ii1), p(ii1), q(ii1)
            end if
         end do
         call draw (p, 1.0_real64, q, cbz, mtp0)
!
      else
!
!   Mode (2)
!
         do k = 1, orbital_count
            write (99, '("1",1i2,1a2," orbital:"/,2(" --------- r --------- ------- P (r) -------",&
               &" ------- Q (r) -------"))') orbital_principal_n(k), orbital_symmetry_labels(k)
            mfk = orbital_last_grid_index(k)
            nb2 = mfk/2
            if (2*nb2 == mfk) then
               nrows = nb2
            else
               nrows = nb2 + 1
            end if
            do ii = 1, nrows
               ii1 = ii
               ii2 = ii1 + nrows
               if (ii2 <= mfk) then
                  write (99, '(1p,6(1x,1d21.14))') &
                     radial_grid(ii1), radial_large_components(ii1,k), radial_small_components(ii1,k), &
                     & radial_grid(ii2), radial_large_components(ii2,&
                     k), radial_small_components(ii2,k)
               else if (ii1 <= mfk) then
                  write (99, '(1p,6(1x,1d21.14))') &
                     radial_grid(ii1), radial_large_components(ii1,k), radial_small_components(ii1,k)
               end if
            end do
            call draw (radial_large_components(:max_radial_points,k), 1.0_real64, &
               & radial_small_components(:max_radial_points,k), cbz, orbital_last_grid_index(k))
!
         end do
!
      end if
!
      return
!
      return
!
      end subroutine prwf
