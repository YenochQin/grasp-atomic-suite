!***********************************************************************
!                                                                      *
      subroutine lagcon(orbital_index, mpi_process_count)
!                                                                      *
!   This  routine  includes  the Lagrange multiplier contribution in   *
!   the 'exchange' term.                                               *
!
!   Parameter nprocs is added so that it works for mpi and serial
!   programs.
!                                                                      *
!                                           Last update: 08 Dec 1992   *
!   Modified by Xinghong He                 Last update: 17 Aug 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: orbital_index_radix
      use atomic_calculation_state, only: accy, max_radial_points, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orbital_lagrange_multipliers, only: ecv, iecc, nec
      use scf_radial_potentials, only: max_radial_points, xp, xq, yp
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: orbital_index
      integer , intent(in) :: mpi_process_count
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, iecck, l1, l2, m, mfm, i
      real(real64) :: eps, wb, wa, wari
!-----------------------------------------------
!
!
      if (nec == 0) return
!
!ww      EPS = ACCY*0.1D 00
      eps = accy*0.01_real64
!
!   Add contributions from off-diagonal parameters to exchange
!
      wb = 1.0_real64/(ucf(orbital_index)*speed_of_light_au)/mpi_process_count
      do k = 1, nec
!
!   Decode index
!
         iecck = iecc(k)
         l1 = iecck/orbital_label_radix
         l2 = iecck - orbital_label_radix*l1
!
         if (orbital_index == l1) then
            m = l2
         else if (orbital_index == l2) then
            m = l1
         else
            cycle
         end if
!
         wa = ecv(k)*wb
         if (abs(wa) < eps) cycle
!
!   ADD CONTRIBUTIONS TO EXCHANGE TERMS
!
         mfm = orbital_last_grid_index(m)
         do i = 1, mfm
            wari = wa*radial_grid(i)
            xp(i) = xp(i) + wari*radial_small_components(i,m)
            xq(i) = xq(i) - wari*radial_large_components(i,m)
         end do
!
      end do
!
      return
      end subroutine lagcon
