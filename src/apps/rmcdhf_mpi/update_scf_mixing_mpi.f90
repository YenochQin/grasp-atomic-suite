! Update the mixing coefficients used by the MPI SCF calculation.
module update_scf_mixing_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: update_scf_mixing_mpi
contains

    ! This routine computes the level weights, the generalized occupa-
    ! tion numbers,  and average energy for  EOL calculations;  this
    ! information and the eigenvectors are then printed out.
    !
    ! Last revision: 24 Dec 1992
    ! Block version by Xinghong He. Last revision: 05 Aug 1998
    !
    ! Modified by G. Gaigalas. 05 Feb 2017
    subroutine update_scf_mixing_mpi(weighted_average_energy)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use debug_output_settings, only: ldbpg
        use atomic_calculation_state, only: requested_level_indices, selected_level_count, requested_level_weights, &
            level_average_weights
        use eigenstate_expansions, only: shifted_eigenvalues, csf_mixing_coefficients
        use orbital_csf_basis, only: orbital_count
        use scf_iteration_state, only: orbital_occupation_weights
        use level_symmetry_metadata,          only: packed_block_symmetries
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
        use symmetry_block_offsets, only: csf_offsets_by_block, mixing_offsets_by_block
        use block_mean_energies, only: mean_energy_by_block
        use transition_block_indices, only: level_block_indices
        use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_process_count, &
            mpi_rank, mpi_sum

        ! Procedure interfaces
        use csf_density_matrix_element_interface, only: csf_density_matrix_element

        ! Module dependencies
        use csf_quantum_queries, only: csf_occupation

        ! Procedure interfaces
        use write_csf_weights_interface, only: write_csf_weights
        implicit none (type, external)

        ! Arguments
        real(real64) , intent(out) :: weighted_average_energy

        ! Local declarations
        integer :: weight_index, orbital_index, symmetry_block_index, csf_index, mixing_coefficient_offset, &
            & selected_level_index, weight_block_index
        real(real64) :: requested_weight, level_energy, weight_sum, occupation_sum
        logical :: use_extended_optimal_level
        ! Compute weighting factors
        weight_sum = 0.0_real64
        do weight_index = 1, selected_level_count
            requested_weight = requested_level_weights(weight_index)
            if (requested_weight == (-2.0_real64)) then
                level_average_weights(weight_index) = 1.0_real64
            else if (requested_weight == (-1.0_real64)) then
                ! WT(J) = IATJPO(J)
                weight_block_index = level_block_indices(weight_index)     ! Block number of this state
                level_average_weights(weight_index) = abs(packed_block_symmetries(weight_block_index))
            else
                level_average_weights(weight_index) = requested_level_weights(weight_index)
            end if
            weight_sum = weight_sum + level_average_weights(weight_index)
        end do

        level_average_weights(:selected_level_count) = level_average_weights(:selected_level_count)/weight_sum
        ! Compute generalised occupation numbers
        ! <----- Distributed ----->
        use_extended_optimal_level = .true.

        do orbital_index = 1, orbital_count
            occupation_sum = 0.0_real64
            do symmetry_block_index = 1, symmetry_block_count
                do csf_index = mpi_rank + 1, csf_counts_by_block(symmetry_block_index), mpi_process_count
                    occupation_sum = occupation_sum + csf_density_matrix_element(use_extended_optimal_level, &
                        & csf_index,csf_index,symmetry_block_index)*csf_occupation(orbital_index, &
                        & csf_index + csf_offsets_by_block(symmetry_block_index))
                end do
            end do
            ! UCF(J) = SUM
            call mpi_allreduce (occupation_sum, orbital_occupation_weights(orbital_index), 1, mpi_double_precision, &
                mpi_sum, mpi_comm_world, mpi_error_code)
        end do
        ! Write out level energies and weights
        if(mpi_rank == 0) write (*, '(/,"Optimise on the following level(s):"/)')
        weighted_average_energy = 0.0_real64
        mixing_coefficient_offset = 0
        do selected_level_index = 1, selected_level_count
            symmetry_block_index = level_block_indices(selected_level_index)                       ! Block number of this state
            level_energy = mean_energy_by_block(symmetry_block_index) + shifted_eigenvalues(selected_level_index)
            if (mpi_rank == 0) then
                write (*, '("Level ",1i2,4x,"Energy = ",1p,1d19.12,4x,"Weight = ",1d12.5)') &
                    & requested_level_indices(selected_level_index), level_energy, &
                    level_average_weights(selected_level_index)
                if (ldbpg(5)) then
                    call ensure_debug_output()
                    write (angular_debug_unit, *) selected_level_index, symmetry_block_index, &
                        & csf_counts_by_block(symmetry_block_index), &
                        & mixing_offsets_by_block(symmetry_block_index)
                    call ensure_debug_output()
                    write (angular_debug_unit, '(/,"Configuration mixing coefficients:")')
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1x,1p,6d12.4)') (csf_mixing_coefficients(csf_index + &
                        & mixing_coefficient_offset),csf_index=1, &
                        & csf_counts_by_block(symmetry_block_index))
                    mixing_coefficient_offset = mixing_coefficient_offset + csf_counts_by_block(symmetry_block_index)
                end if
            end if
            weighted_average_energy = weighted_average_energy + level_average_weights(selected_level_index)*level_energy
        end do
        ! ~~~~~~~~~
        if (mpi_rank == 0) then
            ! ~~~~~~~~~

            call write_csf_weights (.true.)
            ! Write out average energy
            if (selected_level_count > 1) write (*, '(/,"Weighted average energy of these levels = ",1p,d18.10)') &
                & weighted_average_energy
            write (error_unit,"(A25,1PD18.10)") "Weighted average energy: ", weighted_average_energy
            ! Write out generalized occupation numbers
            write (*, '(/,"Generalised occupation numbers:"/)')
            write (*, '(1x,1p,6d12.4)') (orbital_occupation_weights(csf_index),csf_index=1,orbital_count)

            ! ~~~~~~~~~
        end if
        ! ~~~~~~~~~

        return
    end subroutine update_scf_mixing_mpi
end module update_scf_mixing_mpi_interface
