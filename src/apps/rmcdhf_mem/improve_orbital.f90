!***********************************************************************
subroutine improve_orbital(extended_optimal_level, orbital_index, lsort, dampmx)
    ! Module dependencies
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use correlation_orbital_flags, only: lcorre, max_orbitals
    use scf_damping_factors, only: max_orbitals, mixing_damping_factors, orbital_damping_factors
    use atomic_calculation_state, only: accy, max_radial_points, nsic
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use radial_dirac_workspace, only: max_radial_points, mtp0, p, p0, q, q0, tf, tg, xu, xv
    use serial_mpi_context, only: mpi_rank
    use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_energy_parameters, orbital_origin_exponents, &
        orbital_principal_n, orbital_symmetry_labels, ped
    use scf_orthogonality_control, only: orthst
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii

    ! Procedure interfaces
    use build_orbital_potential_interface, only: build_orbital_potential
    use calculate_deferred_corrections_interface, only: calculate_deferred_corrections
    use solve_radial_orbital_interface, only: solve_radial_orbital

    ! Module dependencies
    use radial_orbital_orthogonalization, only: orthsc

    ! Procedure interfaces
    use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
    use update_scf_mixing_interface, only: update_scf_mixing
    use calculate_lagrange_multipliers_interface, only: calculate_lagrange_multipliers

    ! Module dependencies
    use radial_quadrature, only: quad

    ! Procedure interfaces
    use evaluate_radial_self_consistency_interface, only: evaluate_radial_self_consistency
    use check_orbital_damping_interface, only: check_orbital_damping
    use damp_orbital_update_interface, only: damp_orbital_update
    use orthogonalize_same_kappa_orbitals_interface, only: orthogonalize_same_kappa_orbitals
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: orbital_index
    real(real64) , intent(inout) :: dampmx
    logical, intent(in) :: extended_optimal_level
    logical :: lsort

    ! Local constants
    real(real64), parameter :: p2 = 2.0e-1_real64
    real(real64), parameter :: p005 = 5.0e-3_real64
    real(real64), parameter :: p0001 = 1.0e-4_real64

    ! Local declarations
    integer :: ipr, npts, orbital_sign_reversal_count, jp, radial_node_count, i, nwww
    real(real64) :: ed1, gamaj, ed2, wtaev, dnorm, dnfac, del1, del2, orbital_damping_factor
    logical :: fail, first

    ! Initial values
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
    gamaj = orbital_origin_exponents(orbital_index)
    !
    !   C Froese Fischer's parameters IPR, ED1, ED2 are set and
    !   used in this routine and in DAMPCK
    !
    do
        ed2 = orbital_energy_parameters(orbital_index)
        ed1 = ped(orbital_index)
        !
        !   Set up the exchange potential and arrays XU, XV as appropriate
        !
        !   Set coefficients for YPOT, XPOT, DACON
        !   Compute direct potential, exchange potential
        !   Add in Lagrange-multiplier contribution
        !   Add in derivative-terms contribution
        !
        npts = radial_point_count
        call build_orbital_potential (extended_optimal_level, orbital_index, npts)
        !
        !   Calculate deferred corrections
        !
        call calculate_deferred_corrections (orbital_index)
        !
        !   Solve the Dirac equation
        !
        orbital_sign_reversal_count = 0
        call solve_radial_orbital (orbital_index, fail, orbital_sign_reversal_count, jp, radial_node_count)
        !
        !   Upon failure issue message; take corrective action if possible
        !
        if (fail) then
            if (mpi_rank == 0) write (*, '(/," Failure; equation for orbital ",1i2,1a2,&
                &" could not be solved using method ",1i1)') &
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), method(orbital_index)
            if (method(orbital_index) /= 2) then
                method(orbital_index) = 2
                !XHH orthsc does not have any argument
                !    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
                call orthsc
                !CFF        ... avoid rediagonalization
                !            IF (EOL) THEN
                !               CALL MATRIX
                !               CALL NEWCO (WTAEV)
                !            ENDIF
                call calculate_lagrange_multipliers (extended_optimal_level)
                cycle
            else
                if (mpi_rank == 0) write (*, '(/,/," ****** Error in SUBROUTINE IMPROV ******"/,&
                    &" Convergence not obtained"/)')
                !CALL TIMER (0)
                stop
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

        call evaluate_radial_self_consistency (orbital_index)
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
        del1 = abs(1.0_real64 - ed2/orbital_energy_parameters(orbital_index))
        if (method(orbital_index) == 1) then
            del2 = max(abs(1.0_real64 - sqrt(dnorm)),abs(dnfac - 1.0_real64))
            if (del1<p005 .and. del2>p2) then
                method(orbital_index) = 2
                cycle
            end if
        else
            if (del1<p0001 .and. nsic>1) nsic = nsic - 1
        end if
        exit
    end do
    !
    !   Damp the orbital --- if not converged
    !
    if (scnsty(orbital_index) > accy) then
        call check_orbital_damping (ipr, orbital_index, ed1, ed2)
        orbital_damping_factor = abs(orbital_damping_factors(orbital_index))
    else
        orbital_damping_factor = 0.0_real64                           ! take the whole new orbital
    end if
    call damp_orbital_update (orbital_index, orbital_sign_reversal_count, orbital_damping_factor)

    !   Orthogonalize all orbitals of the same kappa in the order
    !   fixed, spectroscopic, correlation orbitals. The order of
    !   orbitals in the latter two classes are sorted according
    !   to their self-consistency and energy.

    if (orthst) then
        !CALL orthor (J, inv)
        nwww = orbital_count
        call orthogonalize_same_kappa_orbitals (nwww, orbital_index, lsort)
    end if
    !
    !   Print details of iteration
    !
    if (mpi_rank == 0)                                                &
        write (*, '(1x,1i2,1a2,1p,1d16.7,1x,1i2,d11.3,1d10.2,1d10.2,0p,f6.3,1x,1i5,1x,1i5,1x,1i2,1x,1i2)') &
        orbital_principal_n(orbital_index),orbital_symmetry_labels(orbital_index),orbital_energy_parameters(orbital_index), &
        & method(orbital_index), &
        & orbital_origin_parameters(orbital_index),scnsty(orbital_index), &
        !cjb DNORM-1 -> SQRT(DNORM)-1
        !cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
        sqrt(dnorm)-1,orbital_damping_factor,jp,orbital_last_grid_index(orbital_index), &
        & orbital_sign_reversal_count,radial_node_count
    dampmx = max(dampmx,abs(orbital_damping_factor))

    return
end subroutine improve_orbital
