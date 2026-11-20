!***********************************************************************
!                                                                      *
      subroutine dampor(j, inv, odampj)
!                                                                      *
!   This subroutine damps the orbital wave function with index J. it   *
!   also  stores  the  previous  determination  of this orbital.       *
!                                                                      *
!   Call(s) to: [LIB92]: COUNT, RINT.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
      use atomic_calculation_state, only: accy, max_radial_points
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
      use orbital_inversion_flags, only: max_orbitals, noinvt
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
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
      integer, intent(in) :: j
      integer, intent(inout) :: inv
      real(real64), intent(in) :: odampj
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtpo, mtpn, mtp, i, nncff, mfj
      real(real64) :: eps, factor, poldi, qoldi, dnorm, dnfac, sgn
      logical :: check
!-----------------------------------------------
!
!
!   Initialization
!
!ww      EPS = 0.1D 00*ACCY
      eps = 0.01_real64*accy
      check = .not.noinvt(j)
!
!   Damp orbital J using the damping factor  ABS (ODAMP(J)):  ODAMP(J)
!   is restricted to the open interval (-1,1) by  DATSCF ; the meaning
!   of a  negative  ODAMP(J)  is that  ABS (ODAMP(J))  is the constant
!   damping factor;  a positive  ODAMP(J)  can be reset by  subroutine
!   dampck
!
!   Store previous determination of this orbital in arrays P, Q
!
!XHH odampj goes to the argument
!      ODAMPJ = ABS (ODAMP(J))
!
      if (odampj > eps) then
!
         factor = 1.0_real64 - odampj
!
         orbital_origin_parameters(j) = factor*p0 + odampj*orbital_origin_parameters(j)
!
         mtpo = orbital_last_grid_index(j)
         mtpn = mtp0
         mtp0 = mtpo
!
         mtp = max(mtpn,mtpo)
         do i = 1, mtp
            poldi = radial_large_components(i,j)
            radial_large_components(i,j) = factor*p(i) + odampj*radial_large_components(i,j)
            p(i) = poldi
            qoldi = radial_small_components(i,j)
            radial_small_components(i,j) = factor*q(i) + odampj*radial_small_components(i,j)
            q(i) = qoldi
         end do
!
!   Compute normalization factor
!
         orbital_last_grid_index(j) = mtp
         dnorm = rint(j,j,0)
         dnfac = 1.0_real64/sqrt(dnorm)
!
!   Determine if inversion is necessary
!
         if (check) then
            call count (radial_large_components(:max_radial_points,j), mtp, nncff, sgn)
            if (sgn < 0.0_real64) then
               inv = inv + 1
               dnfac = -dnfac
            end if
         end if
!
         orbital_origin_parameters(j) = orbital_origin_parameters(j)*dnfac
         radial_large_components(:mtp,j) = dnfac*radial_large_components(:mtp,j)
         radial_small_components(:mtp,j) = dnfac*radial_small_components(:mtp,j)
!
!   Find new MF(J)
!
         mfj = mtp + 1
    3    continue
         mfj = mfj - 1
         if (abs(radial_large_components(mfj,j)) < eps) then
            radial_large_components(mfj,j) = 0.0_real64
            radial_small_components(mfj,j) = 0.0_real64
            goto 3
         else
            orbital_last_grid_index(j) = mfj
         end if
!
      else
!
         orbital_origin_parameters(j) = p0
!
         mtpo = orbital_last_grid_index(j)
         mtpn = mtp0
         mtp0 = mtpo
!
         mtp = max(mtpn,mtpo)
         do i = 1, mtp
            poldi = radial_large_components(i,j)
            radial_large_components(i,j) = p(i)
            p(i) = poldi
            qoldi = radial_small_components(i,j)
            radial_small_components(i,j) = q(i)
            q(i) = qoldi
         end do
!
         orbital_last_grid_index(j) = mtp
!
      end if
!
      return
      end subroutine dampor
