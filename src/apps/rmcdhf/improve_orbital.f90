module improve_orbital_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: improve_orbital
contains

    !***********************************************************************
    subroutine improve_orbital(extended_optimal_level, orbital_index, sort_by_self_consistency, maximum_orbital_damping)
        ! Module dependencies
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use scf_damping_factors, only: orbital_damping_factors
        use atomic_calculation_state, only: working_tolerance, orbital_updates_per_scf_iteration
        use radial_grid_state, only: radial_point_count, radial_grid_derivative
        use radial_dirac_workspace, only: trial_orbital_last_grid_index, trial_large_component, trial_large_origin_parameter, &
            trial_small_component
        use serial_mpi_context, only: mpi_rank
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_origin_exponents, orbital_principal_n, &
            orbital_symmetry_labels, previous_orbital_energy_changes
        use scf_orthogonality_control, only: enforce_orbital_orthogonality
        use scf_iteration_state, only: radial_solution_methods, orbital_self_consistency_errors
        use radial_quadrature_workspace, only: integration_last_index, radial_workspace_a
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters

        ! Procedure interfaces
        use build_orbital_potential_interface, only: build_orbital_potential
        use calculate_deferred_corrections_interface, only: calculate_deferred_corrections
        use solve_radial_orbital_interface, only: solve_radial_orbital

        ! Module dependencies
        use radial_orbital_orthogonalization, only: orthogonalize_radial_orbitals

        ! Procedure interfaces
        use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
        use update_scf_mixing_interface, only: update_scf_mixing
        use calculate_lagrange_multipliers_interface, only: calculate_lagrange_multipliers

        ! Module dependencies
        use radial_quadrature, only: integrate_radial_grid

        ! Procedure interfaces
        use evaluate_radial_self_consistency_interface, only: evaluate_radial_self_consistency
        use check_orbital_damping_interface, only: check_orbital_damping
        use damp_orbital_update_interface, only: damp_orbital_update
        use orthogonalize_same_kappa_orbitals_interface, only: orthogonalize_same_kappa_orbitals
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        real(real64) , intent(inout) :: maximum_orbital_damping
        logical, intent(in) :: extended_optimal_level
        logical :: sort_by_self_consistency

        ! Local constants
        real(real64), parameter :: method_switch_norm_error_threshold = 2.0e-1_real64
        real(real64), parameter :: method_switch_energy_change_threshold = 5.0e-3_real64
        real(real64), parameter :: iteration_reduction_energy_threshold = 1.0e-4_real64

        ! Local declarations
        integer :: ipr, potential_grid_point_count, orbital_sign_reversal_count, matching_grid_index, radial_node_count, nwww
        real(real64) :: previous_relative_energy_change, origin_exponent, old_orbital_energy, radial_norm_squared, &
            normalization_factor, relative_energy_change, normalization_error, orbital_damping_factor
        logical :: solve_failed, first

        ! Initial values
        !-----------------------------------------------
        !
        !   C Froese Fischer's IPR and ED1 parameter
        !
        data ipr/ 0/
        data previous_relative_energy_change/ 0.0_real64/
        data first/ .false./
        !
        !
        !-----------------------------------------------------------------------
        origin_exponent = orbital_origin_exponents(orbital_index)
        !
        !   C Froese Fischer's parameters IPR, ED1, ED2 are set and
        !   used in this routine and in DAMPCK
        !
        do
            old_orbital_energy = orbital_energy_parameters(orbital_index)
            previous_relative_energy_change = previous_orbital_energy_changes(orbital_index)
            !
            !   Set up the exchange potential and arrays XU, XV as appropriate
            !
            !   Set coefficients for YPOT, XPOT, DACON
            !   Compute direct potential, exchange potential
            !   Add in Lagrange-multiplier contribution
            !   Add in derivative-terms contribution
            !
            potential_grid_point_count = radial_point_count
            call build_orbital_potential (extended_optimal_level, orbital_index, potential_grid_point_count)
            !
            !   Calculate deferred corrections
            !
            call calculate_deferred_corrections (orbital_index)
            !
            !   Solve the Dirac equation
            !
            orbital_sign_reversal_count = 0
            call solve_radial_orbital (orbital_index, solve_failed, orbital_sign_reversal_count, matching_grid_index, &
                & radial_node_count)
            !
            !   Upon failure issue message; take corrective action if possible
            !
            if (solve_failed) then
                if (mpi_rank == 0) write (*, '(/," Failure; equation for orbital ",1i2,1a2,&
                    &" could not be solved using method ",1i1)') &
                    orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index), &
                        radial_solution_methods(orbital_index)
                if (radial_solution_methods(orbital_index) /= 2) then
                    radial_solution_methods(orbital_index) = 2
                    !XHH orthsc does not have any argument
                    !    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
                    call orthogonalize_radial_orbitals
                    !CFF        ... avoid rediagonalization
                    !            IF (EOL) THEN
                    !               CALL MATRIX
                    !               CALL NEWCO (WTAEV)
                    !            ENDIF
                    call calculate_lagrange_multipliers (extended_optimal_level)
                    cycle
                else
                    if (mpi_rank == 0) write (*, '(/,/," ****** Error in SUBROUTINE improve_orbital ******"/,&
                        &" Convergence not obtained"/)')
                    !CALL TIMER (0)
                    error stop
                end if
            end if
            !
            !   Compute norm of radial function
            !
            radial_workspace_a(1) = 0.0_real64
            radial_workspace_a(2:trial_orbital_last_grid_index) = (trial_large_component(2:trial_orbital_last_grid_index)**2+ &
                trial_small_component(2:trial_orbital_last_grid_index)**2)*radial_grid_derivative(2:trial_orbital_last_grid_index)
            integration_last_index = trial_orbital_last_grid_index

            call integrate_radial_grid (radial_norm_squared)

            !   Determine self-consistency [multiplied by SQRT(UCF(J))]

            call evaluate_radial_self_consistency (orbital_index)
            !
            !   Normalize
            !
            normalization_factor = 1.0_real64/sqrt(radial_norm_squared)
            trial_large_origin_parameter = trial_large_origin_parameter*normalization_factor
            trial_large_component(:trial_orbital_last_grid_index) = trial_large_component(:trial_orbital_last_grid_index)* &
                normalization_factor
            trial_small_component(:trial_orbital_last_grid_index) = trial_small_component(:trial_orbital_last_grid_index)* &
                normalization_factor
            !
            !   Check if different method should be used or if improvement
            !   count should be reduced
            !
            relative_energy_change = abs(1.0_real64 - old_orbital_energy/orbital_energy_parameters(orbital_index))
            if (radial_solution_methods(orbital_index) == 1) then
                normalization_error = max(abs(1.0_real64 - sqrt(radial_norm_squared)),abs(normalization_factor - 1.0_real64))
                if (relative_energy_change<method_switch_energy_change_threshold .and. &
                    & normalization_error>method_switch_norm_error_threshold) then
                    radial_solution_methods(orbital_index) = 2
                    cycle
                end if
            else
                if (relative_energy_change<iteration_reduction_energy_threshold .and. orbital_updates_per_scf_iteration>1) &
                    orbital_updates_per_scf_iteration = orbital_updates_per_scf_iteration - 1
            end if
            exit
        end do
        !
        !   Damp the orbital --- if not converged
        !
        if (orbital_self_consistency_errors(orbital_index) > working_tolerance) then
            call check_orbital_damping (ipr, orbital_index, previous_relative_energy_change, old_orbital_energy)
            orbital_damping_factor = abs(orbital_damping_factors(orbital_index))
        else
            orbital_damping_factor = 0.0_real64                           ! take the whole new orbital
        end if
        call damp_orbital_update (orbital_index, orbital_sign_reversal_count, orbital_damping_factor)

        !   Orthogonalize all orbitals of the same kappa in the order
        !   fixed, spectroscopic, correlation orbitals. The order of
        !   orbitals in the latter two classes are sorted according
        !   to their self-consistency and energy.

        if (enforce_orbital_orthogonality) then
            !CALL orthor (J, inv)
            nwww = orbital_count
            call orthogonalize_same_kappa_orbitals (nwww, orbital_index, sort_by_self_consistency)
        end if
        !
        !   Print details of iteration
        !
        if (mpi_rank == 0)                                                &
            write (*, '(1x,1i2,1a2,1p,1d16.7,1x,1i2,d11.3,1d10.2,1d10.2,0p,f6.3,1x,1i5,1x,1i5,1x,1i2,1x,1i2)') &
            orbital_principal_n(orbital_index),orbital_symmetry_labels(orbital_index),orbital_energy_parameters(orbital_index), &
            & radial_solution_methods(orbital_index), &
            & orbital_origin_parameters(orbital_index),orbital_self_consistency_errors(orbital_index), &
            !cjb DNORM-1 -> SQRT(DNORM)-1
            !cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
            sqrt(radial_norm_squared)-1,orbital_damping_factor,matching_grid_index,orbital_last_grid_index(orbital_index), &
            & orbital_sign_reversal_count,radial_node_count
        maximum_orbital_damping = max(maximum_orbital_damping,abs(orbital_damping_factor))

        return
    end subroutine improve_orbital
end module improve_orbital_interface
