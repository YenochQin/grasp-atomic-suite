!***********************************************************************
!                                                                      *
    real(real64)  function vpintf (ia,ib)
!                                                                      *
!   Computes nuclear vacuum polarization integrals.                    *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
     use, intrinsic :: iso_fortran_env, only: real64
     use grasp_parameters,   only: orbital_index_radix
     use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
     use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
         radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
     use nuclear_charge_distribution, only: max_radial_points, zdist
     use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
     use radial_quadrature_workspace, only: mtp, ta, tb
     use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
         pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
     use radial_quadrature, only: quad
     implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ia, ib
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       integer :: k
!------------------------------------------------
!
      mtp = min (orbital_last_grid_index(ia),orbital_last_grid_index(ib))
      !TA(1) = 0.0D 00  !! Why would this be wrong?!!
      ta(1)  = 0
      do 1 k = 2,mtp
         ta(k) = (radial_large_components(k,ia)*radial_large_components(k,ib)+radial_small_components(k,ia)* &
            & radial_small_components(k,ib))                  &
                 *zdist(k)
    1 continue
      call quad (vpintf)
!
      if (ldbpr(9)) write (99,'(/"VPINTF: V (",1i2,1a2,",",1i2,1a2,") = ",1pd19.12)') &
         orbital_principal_n(ia),orbital_symmetry_labels(ia),orbital_principal_n(ib), &
         & orbital_symmetry_labels(ib),vpintf
!
      return
!
!
      end function vpintf
