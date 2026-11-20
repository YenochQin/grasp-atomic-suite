!***********************************************************************
      subroutine improv(eol, j, lsort, dampmx)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use correlation_orbital_flags, only: lcorre, max_orbitals
      use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
      use atomic_calculation_state, only: accy, max_radial_points, nsic
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
      use mpi_s, only: mpi_rank
      use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_energy_parameters, orbital_origin_exponents, &
          orbital_principal_n, orbital_symmetry_labels, ped
      use scf_orthogonality_control, only: orthst
      use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
          nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cofpot_i, only: cofpot
      use defcor_i, only: defcor
      use solve_i, only: solve
      use radial_orbital_orthogonalization, only: orthsc
      use matrix_i, only: matrix
      use newco_i, only: newco
      use setlag_i, only: setlag
      use radial_quadrature, only: quad
      use consis_i, only: consis
      use dampck_i, only: dampck
      use dampor_i, only: dampor
      use orthy_i, only: orthy
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      real(real64) , intent(inout) :: dampmx
      logical, intent(in) :: eol
      logical :: lsort
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: p2 = 2.0e-1_real64
      real(real64), parameter :: p005 = 5.0e-3_real64
      real(real64), parameter :: p0001 = 1.0e-4_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipr, npts, inv, jp, nnp, i, nwww
      real(real64) :: ed1, gamaj, ed2, wtaev, dnorm, dnfac, del1, del2, odampj
      logical :: fail, first
!-----------------------------------------------
!
!   C Froese Fischer's IPR and ED1 parameter
!
      data ipr/ 0/
      data ed1/ 0.0_real64/
      data first/ .false./
!
!
!-----------------------------------------------------------------------
      gamaj = orbital_origin_exponents(j)
!
!   C Froese Fischer's parameters IPR, ED1, ED2 are set and
!   used in this routine and in DAMPCK
!
    1 continue
      ed2 = orbital_energy_parameters(j)
      ed1 = ped(j)
!
!   Set up the exchange potential and arrays XU, XV as appropriate
!
!   Set coefficients for YPOT, XPOT, DACON
!   Compute direct potential, exchange potential
!   Add in Lagrange-multiplier contribution
!   Add in derivative-terms contribution
!
      npts = radial_point_count
      call cofpot (eol, j, npts)
!
!   Calculate deferred corrections
!
      call defcor (j)
!
!   Solve the Dirac equation
!
      inv = 0
      call solve (j, fail, inv, jp, nnp)
!
!   Upon failure issue message; take corrective action if possible
!
      if (fail) then
         if (mpi_rank == 0) write (*, '(/," Failure; equation for orbital ",1i2,1a2,&
            &" could not be solved using method ",1i1)') &
            orbital_principal_n(j), orbital_symmetry_labels(j), method(j)
         if (method(j) /= 2) then
            method(j) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            call orthsc
!CFF        ... avoid rediagonalization
!            IF (EOL) THEN
!               CALL MATRIX
!               CALL NEWCO (WTAEV)
!            ENDIF
            call setlag (eol)
            goto 1
         else
            if (mpi_rank == 0) write (*, '(/,/," ****** Error in SUBROUTINE IMPROV ******"/,&
               &" Convergence not obtained"/)')
            !CALL TIMER (0)
            error stop
         end if
      end if
!
!   Compute norm of radial function
!
      ta(1) = 0.0_real64
      ta(2:mtp0) = (p(2:mtp0)**2+q(2:mtp0)**2)*radial_grid_derivative(2:mtp0)
      mtp = mtp0

      call quad (dnorm)

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      call consis (j)
!
!   Normalize
!
      dnfac = 1.0_real64/sqrt(dnorm)
      p0 = p0*dnfac
      p(:mtp0) = p(:mtp0)*dnfac
      q(:mtp0) = q(:mtp0)*dnfac
!
!   Check if different method should be used or if improvement
!   count should be reduced
!
      del1 = abs(1.0_real64 - ed2/orbital_energy_parameters(j))
      if (method(j) == 1) then
         del2 = max(abs(1.0_real64 - sqrt(dnorm)),abs(dnfac - 1.0_real64))
         if (del1<p005 .and. del2>p2) then
            method(j) = 2
            goto 1
         end if
      else
         if (del1<p0001 .and. nsic>1) nsic = nsic - 1
      end if
!
!   Damp the orbital --- if not converged
!
      if (scnsty(j) > accy) then
         call dampck (ipr, j, ed1, ed2)
         odampj = abs(orbital_damping_factors(j))
      else
         odampj = 0.0_real64                           ! take the whole new orbital
      end if
      call dampor (j, inv, odampj)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      if (orthst) then
         !CALL orthor (J, inv)
         nwww = orbital_count
         call orthy (nwww, j, lsort)
      end if
!
!   Print details of iteration
!
      if (mpi_rank == 0)                                                &
         write (*, '(1x,1i2,1a2,1p,1d16.7,1x,1i2,d11.3,1d10.2,1d10.2,0p,f6.3,1x,1i5,1x,1i5,1x,1i2,1x,1i2)') &
            orbital_principal_n(j),orbital_symmetry_labels(j),orbital_energy_parameters(j),method(j), &
            & orbital_origin_parameters(j),scnsty(j), &
!cjb DNORM-1 -> SQRT(DNORM)-1
!cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
                        sqrt(dnorm)-1,odampj,jp,orbital_last_grid_index(j),inv,nnp
      dampmx = max(dampmx,abs(odampj))

      return
      end subroutine improv
