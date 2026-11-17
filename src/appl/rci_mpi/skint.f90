!***********************************************************************
!                                                                      *
      real(real64) function skint (rac, ia, ic, rbd, ib, id, k, iw)
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
      use parameter_def,   only: max_radial_points
      use bess_c, only: bessj, bessn, bj, double, max_radial_points, tc, td, wij
      use debug_c, only: cutoff, double, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use grid_c, only: radial_grid, radial_grid_derivative, radial_derivative_over_radius
      use orb_c, only: double, max_orbitals, orbital_principal_n, orbital_symmetry_labels
      use tatb_c, only: mtp, ta, tb
      use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use zkf_i, only: zkf
      use quad_i, only: quad
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia, ib, ic, id
      integer, intent(in) :: k, iw
      real(real64), intent(in) :: rac(max_radial_points)
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
      wk = dble(k + k + 1)
      mxrac = min(orbital_last_grid_index(ia),orbital_last_grid_index(ic))
      mxrbd = min(orbital_last_grid_index(ib),orbital_last_grid_index(id))
!
!           (k-1)
!  Compute Z     (rho  ; s)
!                    ac
!
      ta(:mxrac) = rac(:mxrac)
      mtp = mxrac
      call zkf (k - 1, ia, ic)
!
      if (dabs(w) < epsi) then
!
!   W = 0 case
!
         tkeep(:mtp) = tb(:mtp)
         call zkf (k + 1, ia, ic)
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
         call zkf (k - 1, ia, ic)
!
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0_real64
         ta(2:mtp) = ((1.0_real64 + bessn(2,iw,2:mtp))*tb(2:mtp)-tkeep(2:mtp)*bessn&
            (2,iw,2:mtp))*rbd(2:mtp)/radial_grid(2:mtp)**2*radial_derivative_over_radius(2:mtp)
         call quad (valu)
         skint = (wk/w)**2*valu
         ta(:mxrbd) = rbd(:mxrbd)*(1.0_real64 + bessj(2,iw,:mxrbd))
         mtp = mxrbd
         call zkf (k + 1, ib, id)
         mtp = min(mtp,mxrac)
         ta(1) = 0.0_real64
         ta(2:mtp) = rac(2:mtp)*(1.0_real64 + bessn(1,iw,2:mtp))*tb(2:mtp)*radial_grid(2:mtp)&
            *radial_grid_derivative(2:mtp)
         call quad (valu)
         skint = skint - valu*w*w/dble((2*k + 3)*(2*k - 1))
!
      end if
!
      if (ldbpr(11)) write (99, 300) k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ic), &
         & orbital_symmetry_labels(ic), orbital_principal_n(ib)&
         , orbital_symmetry_labels(ib), orbital_principal_n(id), orbital_symmetry_labels(id), iw, valu
!
      return
!
  300 format("  (",1i2,")"/,"S     (",1i2,1a2,",",1i2,1a2,"|",1i2,1a2,",",1i2,&
         1a2,";",1i2,") = ",1p,d19.12)
      return
!
      end function skint
