
!***********************************************************************
!                                                                      *
      subroutine dacon
!                                                                      *
!   This  routine  includes  the  contribution from the off-diagonal   *
!   I(a,b) integrals in the 'exchange' term.                           *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT.                                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: max_radial_points, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_potential_state, only: first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
      use orbital_csf_basis, only: max_orbitals, orbital_kappa
      use scf_radial_potentials, only: max_radial_points, xp, xq, yp
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use radial_wavefunction_derivatives, only: dpbdt
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k, iorb, mfi, i
      real(real64) :: twoc, coeff, fk, rporii, pfi, qfi, zbci
!-----------------------------------------------
!
!
      twoc = speed_of_light_au + speed_of_light_au
!
      do k = 1, ndcof
!
         iorb = nda(k)
         call dpbdt (iorb)
         mfi = orbital_last_grid_index(iorb)
         coeff = da(k)
         fk = dble(orbital_kappa(iorb))
!
         do i = 2, mfi
            rporii = 1.0_real64/(transformed_grid_step*radial_derivative_over_radius(i))
            pfi = radial_large_components(i,iorb)
            qfi = radial_small_components(i,iorb)
            zbci = negative_radius_times_potential(i)/speed_of_light_au
            xp(i) = xp(i) + coeff*(ta(i)*rporii+fk*pfi-(twoc*radial_grid(i)+zbci)*qfi)
            xq(i) = xq(i) + coeff*(tb(i)*rporii-fk*qfi+zbci*pfi)
         end do
!
      end do
!
      return
      end subroutine dacon
