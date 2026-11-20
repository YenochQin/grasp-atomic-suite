!***********************************************************************
!                                                                      *
      real(real64) function rkint (rac, ia, ic, rbd, ib, id, k, iw)
!                                                                      *
!   This routine evaluates the transverse interaction integrals.  If   *
!   IW = 0, it calulates the U(r1,r2) integral; otherwise, it calcu-   *
!   lates  R bar (k; a c | b d ; w)  with  w = wac  if  IW = 1,  and   *
!   w = wbd  if IW = 2.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!               [RCI92]: ZKF.                                          *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use bessel_function_workspace, only: bessj, bessn, bj, max_radial_points, tc, td, wij
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use radial_grid_state, only: radial_derivative_over_radius
      use orbital_csf_basis, only: max_orbitals, orbital_principal_n, orbital_symmetry_labels
      use radial_quadrature_workspace,  only: mtp, ta, tb
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
      integer, intent(in) :: ia, ic
      integer, intent(in) :: ib, id
      integer, intent(in) :: k
      integer, intent(in) :: iw
      real(real64), intent(in) :: rac(max_radial_points)
      real(real64), intent(in) :: rbd(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mxrbd, mxrac, i
      real(real64) :: result
!-----------------------------------------------
!
      mxrbd = min (orbital_last_grid_index(ib),orbital_last_grid_index(id))
      mxrac = min (orbital_last_grid_index(ia),orbital_last_grid_index(ic))
!
      if (iw == 0) then
!
!   IW = 0
!
         ta(:mxrac) = rac(:mxrac)
         mtp = mxrac
         call zkf (k, ia, ic)
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0_real64
         ta(2:mtp) = rbd(2:mtp)*tb(2:mtp)*radial_derivative_over_radius(2:mtp)
!
      else
!
!   IW = 1,2
!
         ta(:mxrac) = rac(:mxrac)*(1.0_real64 + bessj(1,iw,:mxrac))
         mtp = mxrac
         call zkf (k, ia, ic)
         mtp = min(mtp,mxrbd)
         ta(1) = 0.0_real64
         ta(2:mtp) = rbd(2:mtp)*(1.0_real64 + bessn(1,iw,2:mtp))*tb(2:mtp)*radial_derivative_over_radius(2:&
            mtp)
!
      end if
!
      call quad (result)
      rkint = result
!
!   Debug printout if option set
!
      if (ldbpr(11)) write (99, '("_ (",1i2,")"/,"R     (",1i2,1a2,",",1i2,1a2,"|",1i2,1a2,",",1i2,1a2,";",1i2,&
         &") = ",1p,d19.12)') k, orbital_principal_n(ia), orbital_symmetry_labels(ia), orbital_principal_n(ic), &
         & orbital_symmetry_labels(ic), orbital_principal_n(ib)&
         , orbital_symmetry_labels(ib), orbital_principal_n(id), orbital_symmetry_labels(id), iw, result
!
      return
!
      return
!
      end function rkint
