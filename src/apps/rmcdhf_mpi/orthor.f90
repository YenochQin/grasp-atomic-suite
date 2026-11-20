!***********************************************************************
!                                                                      *
      subroutine orthor(orbital_index, orbital_sign_reversal_count)
!
!                                                                      *
!   This routine Schmidt orthogonalizes orbital  J  to  all orbitals   *
!   which  have  better  self-consistency.  Note that fixed orbitals   *
!   have the best self-consistency.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
! INV not initialized here but used as
! INV = INV + 1
! A bug regarding IORDER fixed and new schemes introduced. See comments
! Anyway this routine was not used anywhere.
!XHH 1997.02.21
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      use atomic_calculation_state, only: accy, max_radial_points
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use orbital_overlap_control, only: novl
      use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa
      use scf_orbital_order, only: iorder, max_orbitals
      use orbital_inversion_flags, only: max_orbitals, noinvt
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use weighted_radial_overlap, only: rint
!      USE count_I
      use radial_node_count, only: count
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: orbital_index
      integer , intent(inout) :: orbital_sign_reversal_count
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kraw, k, mtp, i, nncff, jgg
      real(real64) :: eps, ovrlap, original_norm_squared, factor, sgn
      logical :: check, changed
!-----------------------------------------------
!
!XHH Added common /orba/ and /iounit/

!
!ww      EPS = ACCY*0.1D 00
      eps = accy*0.01_real64
!
      check = .not.noinvt(orbital_index)
!
!XHH
! Bug fixed. Here J is actually IORDER(J_raw), thus K should be
! treated the same way.
      changed = .false.
!      DO 4 K = 1,NW
      do kraw = 1, orbital_count
         k = iorder(kraw)
!         write(istde,*) '***',kraw,k,np(k),nh(k),scnsty(k),'***'
!
!XHH orbitals with higher self-consistency are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) .AND.
!     :        (SCNSTY(K) .LT. SCNSTY(J)) ) THEN
!XHH All orbitals are considered
!         IF ( (K .NE. J) .AND.
!     :        (NAK(K) .EQ. NAK(J)) ) THEN
!XHH orbitals before the current and unchanged ones are considered
         if (.not.(orbital_kappa(k)==orbital_kappa(orbital_index) .and. (k<orbital_index .or. is_orbital_fixed(k)))) cycle

         changed = .true.
!
!   Compute overlap
!
         ovrlap = rint(orbital_index,k,0)
!
!   Schmidt orthogonalise
!
         orbital_origin_parameters(orbital_index) = orbital_origin_parameters(orbital_index) - ovrlap*orbital_origin_parameters(k)
         mtp = max(orbital_last_grid_index(orbital_index),orbital_last_grid_index(k))
         radial_large_components(:mtp,orbital_index) = radial_large_components(:mtp, &
            & orbital_index) - ovrlap*radial_large_components(:mtp,k)
         radial_small_components(:mtp,orbital_index) = radial_small_components(:mtp, &
            & orbital_index) - ovrlap*radial_small_components(:mtp,k)
      end do
!XHH
      if (changed) then
!
!   Normalise
!
         orbital_last_grid_index(orbital_index) = mtp
         original_norm_squared = rint(orbital_index,orbital_index,0)
         factor = 1.0_real64/sqrt(original_norm_squared)
!
!   Determine if inversion is necessary
!
         if (check) then
            call count (radial_large_components(:max_radial_points,orbital_index), mtp, nncff, sgn)
            if (sgn < 0.0_real64) then
               orbital_sign_reversal_count = orbital_sign_reversal_count + 1
               factor = -factor
            end if
         end if
!
!   Perform normalization and/or inversion
!
         orbital_origin_parameters(orbital_index) = factor*orbital_origin_parameters(orbital_index)
         radial_large_components(2:mtp,orbital_index) = factor*radial_large_components(2:mtp,orbital_index)
         radial_small_components(2:mtp,orbital_index) = factor*radial_small_components(2:mtp,orbital_index)
!
!   Find new MF(J)
!
         mtp = mtp + 1
    3    continue
         mtp = mtp - 1
         if (abs(radial_large_components(mtp,orbital_index)) < eps) then
            radial_large_components(mtp,orbital_index) = 0.0_real64
            radial_small_components(mtp,orbital_index) = 0.0_real64
            goto 3
         else
            orbital_last_grid_index(orbital_index) = mtp
         end if
!
!XHH two lines moved forward and added an ENDIF
!         ENDIF
!    4 CONTINUE
      end if
!
      return
      end subroutine orthor
