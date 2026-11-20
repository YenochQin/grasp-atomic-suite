!***********************************************************************
!                                                                      *
      real(real64) function skint (recoupling_coefficient, orbital_a, orbital_c, rbd, orbital_b, orbital_d, multipole_rank, iw)
!                                                                      *
!   This routine evaluates transverse interaction integrals:           *
!                                                                      *
!                          (k)                                         *
!                         S   (a,c;b,d;w)                              *
!                                                                      *
!   where w = wac if IW = 1, and w= wbd if IW = 2.                     *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: ZKF.                                          *
!                                           Last update: 06 Nov 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use bessel_function_workspace, only: bessj, bessn, bj, max_radial_points, tc, td, wij
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use radial_grid_state, only: radial_grid, radial_grid_derivative, radial_derivative_over_radius
      use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
      use radial_quadrature_workspace, only: mtp, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zkf_i, only: zkf
      use radial_quadrature, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d
      integer, intent(in) :: multipole_rank, iw
      real(real64), intent(in) :: recoupling_coefficient(max_radial_points)
      real(real64), intent(in) :: rbd(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxrac, mxrbd, i
      real(real64) :: tkeep(max_radial_points)
      real(real64) :: epsi, w, wk, valu
!-----------------------------------------------
!
!
      epsi = 1.0e-10_real64
!
      w = wij(iw)
      wk = dble(multipole_rank + multipole_rank + 1)
      mxrac = min(orbital_last_grid_index(orbital_a),orbital_last_grid_index(orbital_c))
      mxrbd = min(orbital_last_grid_index(orbital_b),orbital_last_grid_index(orbital_d))
!
!           (k-1)
!  Compute Z     (rho  ; s)
!                    ac
!
      ta(:mxrac) = recoupling_coefficient(:mxrac)
      mtp = mxrac
      call zkf (multipole_rank - 1, orbital_a, orbital_c)
!
      if (abs(w) < epsi) then
!
!   W = 0 case
!
         tkeep(:mtp) = tb(:mtp)
         call zkf (multipole_rank + 1, orbital_a, orbital_c)
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0_real64
         ta(2:mtp) = rbd(2:mtp)*radial_derivative_over_radius(2:mtp)*(tb(2:mtp)-tkeep(2:mtp))
         call quad (valu)
         skint = wk*valu*0.5_real64
!
      else
!
!   Finite w: see I P Grant and B J McKenzie, J Phys B: At Mol Phys,
!   13 (1980) 2671-2681
!
         tkeep(:mtp) = tb(:mtp)
         ta(:mtp) = -ta(:mtp)*bessj(1,iw,:mtp)
         call zkf (multipole_rank - 1, orbital_a, orbital_c)
!
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0_real64
         ta(2:mtp) = ((1.0_real64 + bessn(2,iw,2:mtp))*tb(2:mtp)-tkeep(2:mtp)*bessn&
            (2,iw,2:mtp))*rbd(2:mtp)/radial_grid(2:mtp)**2*radial_derivative_over_radius(2:mtp)
         call quad (valu)
         skint = (wk/w)**2*valu
         ta(:mxrbd) = rbd(:mxrbd)*(1.0_real64 + bessj(2,iw,:mxrbd))
         mtp = mxrbd
         call zkf (multipole_rank + 1, orbital_b, orbital_d)
         mtp = min(mtp,mxrac)
         ta(1) = 0.0_real64
         ta(2:mtp) = recoupling_coefficient(2:mtp)*(1.0_real64 + bessn(1,iw,2:mtp))*tb(2:mtp)*radial_grid(2:mtp)&
            *radial_grid_derivative(2:mtp)
         call quad (valu)
         skint = skint - valu*w*w/dble((2*multipole_rank + 3)*(2*multipole_rank - 1))
!
      end if
!
      if (ldbpr(11)) write (99, '("  (",1i2,")"/,"S     (",1i2,1a2,",",1i2,1a2,"|",1i2,1a2,",",1i2,1a2,";",1i2,&
         &") = ",1p,d19.12)') multipole_rank, orbital_principal_n(orbital_a), orbital_symmetry_labels(orbital_a), &
            & orbital_principal_n(orbital_c), &
         & orbital_symmetry_labels(orbital_c), orbital_principal_n(orbital_b)&
         , orbital_symmetry_labels(orbital_b), orbital_principal_n(orbital_d), orbital_symmetry_labels(orbital_d), iw, valu
!
      return
!
      return
!
      end function skint
