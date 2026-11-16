
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
      use def_c
      use grid_c
      use npot_c, only: double, first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
      use orb_c
      use pote_c, only: double, max_radial_points, xp, xq, yp
      use scf_c
      use tatb_c, only: double, mtp, radial_workspace_capacity, ta, tb
      use wave_c
      use dpbdt_i, only: dpbdt
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
