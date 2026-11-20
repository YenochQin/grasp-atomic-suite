!***********************************************************************
!                                                                      *
      subroutine xpot(orbital_index)
      use debug_output_settings, only: angular_debug_unit, ensure_debug_output
!                                                                      *
!   This subroutine tabulates the exchange terms (the first terms on   *
!   the  right-hand  sides of eqs (14), I P Grant, B J McKenzie, P H   *
!   Norrington, D F Mayers, and N C Pyper, Computer  Phys  Commun 21   *
!   (1980) 211) for orbital J. The exchange terms are stored  in the   *
!   common arrays XP and XQ.                                           *
!                                                                      *
!   Call(s) to: [LIB92]: DRAW, YZK.                                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: orbital_index_radix
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: max_radial_points, speed_of_light_au
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
      use scf_radial_potentials, only: max_radial_points, xp, xq, yp
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use hartree_yz_functions, only: yzk
      use radial_ascii_plot, only: draw
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: orbital_index
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: orbital_label_radix = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: radial_index, index, label, k, ioy1, ioy2, iorb, nb2, row_count, ii, ii1, ii2
      real(real64) :: coeff, ctb
!-----------------------------------------------
!
!   Debug printout: header
!
      if (ldbpr(27) .or. ldbpr(28)) then
         call ensure_debug_output()
         write (angular_debug_unit, '(/,/," Exchange potential contributions (coefficients will ",&
         &" be divided by C) for ",1i2,1a2," orbital :"/,/)') orbital_principal_n(orbital_index), &
            & orbital_symmetry_labels(orbital_index)
      end if
!
!   Clear for accumulation of sums
!
      xp(:radial_point_count) = 0.0_real64
      xq(:radial_point_count) = 0.0_real64
!
!   Add contributions from exchange terms
!
      do index = 1, nxcof

         ! Decode information in label
         label = nxa(index)
         k = mod(label,orbital_label_radix)
         label = label/orbital_label_radix
         ioy1 = mod(label,orbital_label_radix)
         label = label/orbital_label_radix
         ioy2 = mod(label,orbital_label_radix)
         iorb = label/orbital_label_radix
         coeff = xa(index)

         ! Debug printout: composition
         if (ldbpr(27)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,25x,"(",1i2,")"/,1x,1p,d21.14,"* Y    (",1i2,1a2,",",1i2,1a2,") ","* P (",&
            &1i2,1a2,")")') k, coeff, orbital_principal_n(ioy1), orbital_symmetry_labels(ioy1), &
            & orbital_principal_n(ioy2)&
            , orbital_symmetry_labels(ioy2), orbital_principal_n(iorb), orbital_symmetry_labels(iorb)
         end if

         call yzk (k, ioy1, ioy2)
!
!   Accumulate contributions
!
         coeff = coeff/speed_of_light_au
         !DO I = 1, MF(IORB)
         do radial_index = 1, radial_point_count
            ctb = coeff*tb(radial_index)
            xp(radial_index) = xp(radial_index) + ctb*radial_small_components(radial_index,iorb)
            xq(radial_index) = xq(radial_index) - ctb*radial_large_components(radial_index,iorb)
         end do
      end do
!
!   Debug printout: potential functions
!
      if (ldbpr(28)) then
         call ensure_debug_output()
         write (angular_debug_unit, '(/,/,31x,"(P)",19x,"(Q)",41x,"(P)",19x,"(Q)"/,&
            &2(" --------- r --------- ------ X  (r) -------"," ------ X  (r) -------"))')
         nb2 = radial_point_count/2
         if (2*nb2 == radial_point_count) then
            row_count = nb2
         else
            row_count = nb2 + 1
         end if
         do ii = 1, row_count
            ii1 = ii
            ii2 = ii1 + row_count
            if (ii2 <= radial_point_count) then
               call ensure_debug_output()
               write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), xp(ii1), xq(ii1), radial_grid(ii2), xp(ii2), xq(&
                  ii2)
            else if (ii1 <= radial_point_count) then
               call ensure_debug_output()
               write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), xp(ii1), xq(ii1)
            end if
         end do
         call draw (xp, 1.0_real64, xq, speed_of_light_au, radial_point_count)
      end if
!
      return
!
      return
!
      end subroutine xpot
