! Iterate the Hamiltonian solve and radial orbital updates until the SCF stopping criteria are met.
module run_scf_iterations_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: run_scf_iterations
contains

    ! This  subroutine  performs  the SCF iterations. The procedure is
    ! essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3
    ! (1986) 290.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 22 Dec 1992
    ! Block version by Xinghong He. Last revision: 05 Aug 1998
    ! Modified by G. Gaigalas. 05 Feb 2017
    subroutine run_scf_iterations(extended_optimal_level, output_wavefunction_path)
        ! Module dependencies
        use calculation_file_units, only: mixing_output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc, dalloc
        use transition_block_indices, only: level_block_indices
        use calculation_input_defaults, only: nondefault_options_flag
        use atomic_calculation_state, only: working_tolerance, requested_level_indices, selected_level_count, &
            maximum_scf_iterations, orbital_updates_per_scf_iteration, level_average_weights
        use debug_output_settings, only: rmcdhf_debug_flags
        use eigenstate_expansions, only: shifted_eigenvalues, csf_mixing_coefficients
        use fixed_orbital_flags, only: is_orbital_fixed
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use standard_io_units, only: stderr_unit, stdin_unit
        use orbital_lagrange_multipliers, only: lagrange_multiplier_values, lagrange_orbital_pair_labels, lagrange_multiplier_count
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use serial_mpi_context, only: mpi_rank
        use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, csf_offsets_by_block, &
            eigenvalue_offsets_by_block, mixing_offsets_by_block
        use block_mean_energies, only: mean_energy_by_block
        use orbital_csf_basis, only: csf_count, csf_subshell_occupations, orbital_count
        use scf_orbital_order, only: orbital_update_order
        use scf_iteration_state, only: orbital_self_consistency_errors, direct_coefficient_capacity, &
            exchange_coefficient_capacity, y_potential_coefficient_capacity, direct_coefficients, exchange_coefficients, &
            y_potential_coefficients, direct_orbital_indices, exchange_coefficient_labels, y_potential_coefficient_labels
        use scf_orthogonality_control, only: enforce_orbital_orthogonality

        ! Procedure interfaces
        use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
        use update_scf_mixing_interface, only: update_scf_mixing
        use calculate_lagrange_multipliers_interface, only: calculate_lagrange_multipliers
        use improve_orbital_interface, only: improve_orbital
        use find_least_consistent_orbital_interface, only: find_least_consistent_orbital
        use print_radial_wavefunction_interface, only: print_radial_wavefunction

        ! Module dependencies
        use radial_orbital_orthogonalization, only: orthogonalize_radial_orbitals

        ! Procedure interfaces
        use write_scf_radial_wavefunctions_interface, only: write_scf_radial_wavefunctions
        use finish_summary_interface, only: finish_summary
        implicit none (type, external)

        ! Arguments
        !
        !   G l o b a l   P a r a m e t e r s
        logical, intent(in) :: extended_optimal_level
        character  :: output_wavefunction_path*(*)

        ! Local declarations
        integer :: orbital_index, symmetry_block_index, mcp_unit, scf_iteration, orbital_update_index, extra_orbital_update, &
            & least_consistent_orbital
        real(real64) :: weighted_average_energy, previous_weighted_average_energy, maximum_orbital_damping
        logical :: scf_converged, sort_by_self_consistency, first_davidson_solve

        total_csf_count = csf_count

        ! Determine Orthonomalization order --- lsort

        if (nondefault_options_flag == 0) then
            sort_by_self_consistency = .false.
        else
            do
                write (stderr_unit, *) "Orthonomalization order? "
                write (stderr_unit, *) "     1 -- Update order"
                write (stderr_unit, *) "     2 -- Self consistency connected"
                read (stdin_unit, *) orbital_index
                if (orbital_index == 1) then
                    sort_by_self_consistency = .false.
                else if (orbital_index == 2) then
                    sort_by_self_consistency = .true.
                else
                    write (stderr_unit, *) "Input is wrong, redo..."
                    cycle
                end if
                exit
            end do
        end if

        ! Deallocate storage that will no longer be used

        ! Allocate and fill in auxiliary arrays

        call alloc (csf_offsets_by_block, symmetry_block_count, "csf_offsets_by_block", "run_scf_iterations")
        call alloc (eigenvalue_offsets_by_block, symmetry_block_count, "eigenvalue_offsets_by_block", "run_scf_iterations")
        call alloc (mixing_offsets_by_block, symmetry_block_count, "mixing_offsets_by_block", "run_scf_iterations")
        call alloc (mean_energy_by_block, symmetry_block_count, "mean_energy_by_block", "run_scf_iterations")

        csf_offsets_by_block(1) = 0
        eigenvalue_offsets_by_block(1) = 0
        mixing_offsets_by_block(1) = 0
        do symmetry_block_index = 2, symmetry_block_count
            csf_offsets_by_block(symmetry_block_index) = csf_offsets_by_block(symmetry_block_index-1) + &
                & csf_counts_by_block(symmetry_block_index - 1)
            eigenvalue_offsets_by_block(symmetry_block_index) = eigenvalue_offsets_by_block(symmetry_block_index-1) + &
                & eigenvector_counts_by_block(symmetry_block_index - 1)
            mixing_offsets_by_block(symmetry_block_index) = mixing_offsets_by_block(symmetry_block_index-1) + &
                & eigenvector_counts_by_block(symmetry_block_index - 1) &
                & *csf_counts_by_block(symmetry_block_index - 1)
        end do

        ! Size of the eigenvector array for all blocks
        mixing_coefficient_storage_size = mixing_offsets_by_block(symmetry_block_count) + &
            & eigenvector_counts_by_block(symmetry_block_count)*csf_counts_by_block(symmetry_block_count)

        if (extended_optimal_level) then
            call alloc (shifted_eigenvalues, selected_level_count, "shifted_eigenvalues", "run_scf_iterations")
            call alloc (csf_mixing_coefficients, mixing_coefficient_storage_size, "csf_mixing_coefficients", "run_scf_iterations")

        end if

        direct_coefficient_capacity = 0
        exchange_coefficient_capacity = 0
        y_potential_coefficient_capacity = 0

        ! This call should only be made AFTER the call to newco

        ! For (E)OL calculations, determine the level energies and
        ! mixing coefficients

        ! .. set the logical variable dvdfirst
        first_davidson_solve = .true.
        if (extended_optimal_level) then
            call solve_hamiltonian_blocks (first_davidson_solve)
            call update_scf_mixing (weighted_average_energy)
        end if
        previous_weighted_average_energy = 0.0
        first_davidson_solve = .false.
        do scf_iteration = 1, maximum_scf_iterations
            if (mpi_rank == 0) write (*, '(/," Iteration number ",1i3,/," --------------------")') scf_iteration

            ! For all pairs constrained through a Lagrange multiplier, compute
            ! the Lagrange multiplier

            call calculate_lagrange_multipliers (extended_optimal_level)

            ! Improve all orbitals in turn

            maximum_orbital_damping = 0.0
            if (mpi_rank == 0) write (*, '(41x,"Self-            Damping"/,"Subshell    Energy    Method   P0    ",&
                &"consistency  Norm-1  factor  JP"," MTP INV NNP"/)')
            do orbital_index = 1, orbital_count
                orbital_update_index = orbital_update_order(orbital_index)
                if (is_orbital_fixed(orbital_update_index)) cycle
                call improve_orbital (extended_optimal_level, orbital_update_index, sort_by_self_consistency, &
                    & maximum_orbital_damping)
            end do
            ! For KOUNT = 1 to NSIC: find the least self-consistent orbital;
            ! improve it
            do extra_orbital_update = 1, orbital_updates_per_scf_iteration
                call find_least_consistent_orbital (least_consistent_orbital)
                if (least_consistent_orbital == 0) then
                    scf_converged = .true.
                    exit
                else
                    if (orbital_self_consistency_errors(least_consistent_orbital) <= working_tolerance) then
                        scf_converged = .true.
                        exit
                    end if
                end if
                call improve_orbital (extended_optimal_level, least_consistent_orbital, sort_by_self_consistency, &
                    & maximum_orbital_damping)
            end do

            if (extra_orbital_update > orbital_updates_per_scf_iteration) then
                call find_least_consistent_orbital (least_consistent_orbital)

                if (least_consistent_orbital == 0) then
                    scf_converged = .true.
                else
                    if (orbital_self_consistency_errors(least_consistent_orbital) <= working_tolerance) then
                        scf_converged = .true.
                    else
                        scf_converged = .false.
                    end if
                end if

            end if
            if (rmcdhf_debug_flags(24) .and. mpi_rank==0) call print_radial_wavefunction (0)

            ! Perform Gram-Schmidt process
            ! For OL calculation, orthst is true and orbitals are orthonormalized
            ! in subroutine improv. For AL calculation, orthst is false.
            if (.not.enforce_orbital_orthogonality) call orthogonalize_radial_orbitals

            ! Write this iteration's radial wavefunctions to the selected output path.

            if (mpi_rank == 0) call write_scf_radial_wavefunctions (output_wavefunction_path)

            if (extended_optimal_level) then
                call solve_hamiltonian_blocks(first_davidson_solve)
                call update_scf_mixing (weighted_average_energy)
            end if
            !         Make this a relative convergence test
            ! cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
            ! DAMPMX.LT.1.0D-2) CONVG=.true.
            if (abs((weighted_average_energy - previous_weighted_average_energy)/weighted_average_energy) < &
                & 0.001_real64*working_tolerance) scf_converged = .true.
            previous_weighted_average_energy = weighted_average_energy
            if (.not.scf_converged) cycle
            if (rmcdhf_debug_flags(25) .and. .not.rmcdhf_debug_flags(24) .and. mpi_rank==0) call print_radial_wavefunction (0)
            exit

        end do

        if (scf_iteration > maximum_scf_iterations) then
            if (mpi_rank == 0) write (stderr_unit, *) " Maximum iterations in SCF Exceeded."
        end if

        do mcp_unit = 31, 32 + max_mcp_file_multipole_rank
            close(mcp_unit)                                ! The MCP coefficient files
        end do

        if (mpi_rank == 0) then
            close(mixing_output_unit)                               ! The .mix file
            ! Complete the summary - moved from rscf92 for easier alloc/dalloc
            call finish_summary
        end if
        ! Deallocate storage
        call dalloc (level_average_weights, "level_average_weights", "run_scf_iterations")                       !Either getold or getald

        if (lagrange_multiplier_count > 0) then
            call dalloc (lagrange_orbital_pair_labels, "lagrange_orbital_pair_labels", "run_scf_iterations")
            call dalloc (lagrange_multiplier_values, "lagrange_multiplier_values", "run_scf_iterations")
            call dalloc (csf_subshell_occupations, "csf_subshell_occupations", "run_scf_iterations")
        end if

        if (direct_coefficient_capacity > 0) then
            call dalloc (direct_coefficients, "direct_coefficients", "run_scf_iterations")
            call dalloc (direct_orbital_indices, "direct_orbital_indices", "run_scf_iterations")
        end if

        if (exchange_coefficient_capacity > 0) then
            call dalloc (exchange_coefficients, "exchange_coefficients", "run_scf_iterations")
            call dalloc (exchange_coefficient_labels, "exchange_coefficient_labels", "run_scf_iterations")
        end if

        if (y_potential_coefficient_capacity > 0) then
            call dalloc (y_potential_coefficients, "y_potential_coefficients", "run_scf_iterations")
            call dalloc (y_potential_coefficient_labels, "y_potential_coefficient_labels", "run_scf_iterations")
        end if

        if (extended_optimal_level) then
            call dalloc (shifted_eigenvalues, "shifted_eigenvalues", "run_scf_iterations")
            call dalloc (csf_mixing_coefficients, "csf_mixing_coefficients", "run_scf_iterations")

            call dalloc (highest_selected_level_by_block, "highest_selected_level_by_block", "run_scf_iterations") ! getold.f
            call dalloc (mean_energy_by_block, "mean_energy_by_block", "run_scf_iterations")     ! getold.f
            call dalloc (level_block_indices, "level_block_indices", "run_scf_iterations")     ! Allocated in getold.f
            call dalloc (requested_level_indices, "requested_level_indices", "run_scf_iterations")       ! Allocated in items.f<-getold.f
        end if
        call dalloc (csf_offsets_by_block, "csf_offsets_by_block", "run_scf_iterations")
        call dalloc (eigenvalue_offsets_by_block, "eigenvalue_offsets_by_block", "run_scf_iterations")
        call dalloc (mixing_offsets_by_block, "mixing_offsets_by_block", "run_scf_iterations")

        return
    end subroutine run_scf_iterations
end module run_scf_iterations_interface
