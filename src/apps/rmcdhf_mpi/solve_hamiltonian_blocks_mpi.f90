! Distribute Hamiltonian assembly and solve the requested symmetry blocks.
module solve_hamiltonian_blocks_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_hamiltonian_blocks_mpi
contains

    !   Calls routines to form the Hamiltonian matrix and to diagonalise
    !   it. The total angular momenta and parity of each  ASF  is found;
    !   the eigenvectors are normalised so  that the sign of the largest
    !   element is positive.
    !
    ! Last revision: 24 Dec 1992
    ! Block version by Xinghong He. 07 Aug 1998
    subroutine solve_hamiltonian_blocks_mpi(first_davidson_solve)
        ! Module dependencies
        use mpi_mcp_calculation_input, only: indexed_input
        use rmcdhf_mcp_io, only: read_sparse, profile_indexed_matrix_skip, mcp_now, profile_stage, verification_enabled, log_bits
        use mcp_density_contractions, only: record_solver_workspace
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: mixing_output_unit
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use array_memory_management, only: alloc, dalloc
        use scf_damping_factors, only: mixing_damping_factors
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use debug_output_settings, only: ldbpg
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, symmetry_block_count, &
            highest_selected_level_by_block
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
            matrix_partner_indices
        use mpi_state, only: mpi_allreduce, mpi_barrier, mpi_comm_world, mpi_double_precision, mpi_error_code, &
            mpi_process_count, mpi_rank, mpi_sum
        use orbital_csf_basis, only: csf_count, orbital_count
        use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, csf_offsets_by_block, &
            eigenvalue_offsets_by_block, mixing_offsets_by_block
        use block_mean_energies, only: mean_energy_by_block
        use level_symmetry_metadata,          only: packed_block_symmetries

        ! Procedure interfaces
        use assemble_hamiltonian_block_interface, only: assemble_hamiltonian_block
        use diagonalize_hamiltonian_block_mpi_interface, only: diagonalize_hamiltonian_block_mpi
        use blas_lapack_interfaces, only: dcopy, ddot, dscal
        implicit none (type, external)

        ! Arguments
        logical, intent(in) ::  first_davidson_solve

        ! Local declarations
        integer :: file_electron_count, discarded_header_value, symmetry_block_index, i, &
            block_csf_offset, block_eigenvalue_offset, block_mixing_offset, csf_index, idiag, level_index, iofset, &
            jother, ia, jofset, iattmp, iastmp
        integer(int64) :: matrix_entry
        real(real64), pointer :: cmvl(:)
        real(real64) :: amax, cdampj, dnfac, evecij, sum, tmp, omcdaj, ovrlap, wa, started
        logical :: first

        save first

        ! Initial values
        data first/ .true./

        if (mpi_rank == 0) write (output_unit, *)

        ! Allocate memory for CMVL once (the maximum size)
        ! Save previous estimate of eigenvectors

        if (.not.first) then
            call alloc (cmvl, mixing_coefficient_storage_size, "CMVL", "solve_hamiltonian_blocks_mpi")
            call dcopy (mixing_coefficient_storage_size, csf_mixing_coefficients, 1, cmvl, 1)
        end if

        ! Position the files - MCP files (unit NFILE) for reading
        ! and mixing coefficients file (unit 25) for writing

        ! To put in ncmin and nvecsiz. Values read here are the same
        ! as the values already loaded in the shared modules.
        if (mpi_rank == 0) then
            rewind (mixing_output_unit)
            read (mixing_output_unit)                               ! 'G92MIX'
            read (mixing_output_unit) file_electron_count, total_csf_count, orbital_count, discarded_header_value, &
                & discarded_header_value, symmetry_block_count
            backspace (mixing_output_unit)
            write (mixing_output_unit) file_electron_count, total_csf_count, orbital_count, selected_level_count, &
                mixing_coefficient_storage_size, &
                symmetry_block_count
        end if

        ! Do the job block by block

        do symmetry_block_index = 1, symmetry_block_count                      ! block do-loop

            started = mcp_now()
            csf_count = csf_counts_by_block(symmetry_block_index)
            if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
                call profile_indexed_matrix_skip(symmetry_block_index)
                call profile_stage("matrix_skipped",0,symmetry_block_index,started)
                cycle
            end if

            stored_matrix_element_count = indexed_input%row_counts(symmetry_block_index)
            call alloc(matrix_partner_indices,stored_matrix_element_count,"matrix_partner_indices","solve_hamiltonian_blocks_mpi")
            call alloc(stored_segment_end_offsets,0,csf_count,"stored_segment_end_offsets","solve_hamiltonian_blocks_mpi")
            stored_segment_end_offsets = 0
            call read_sparse(symmetry_block_index, &
                stored_segment_end_offsets(mpi_rank+1:csf_count:mpi_process_count),matrix_partner_indices,"matrix")

            block_csf_offset = csf_offsets_by_block(symmetry_block_index)
            block_eigenvalue_offset = eigenvalue_offsets_by_block(symmetry_block_index)
            block_mixing_offset = mixing_offsets_by_block(symmetry_block_index)

            call alloc(hamiltonian_values,stored_matrix_element_count,"hamiltonian_values","solve_hamiltonian_blocks_mpi")
            hamiltonian_values = 0.0_real64
            call record_solver_workspace(csf_count,stored_matrix_element_count, &
                highest_selected_level_by_block(symmetry_block_index),mixing_coefficient_storage_size)
            call assemble_hamiltonian_block (symmetry_block_index, mpi_rank, mpi_process_count)
            if (verification_enabled) then
                do matrix_entry = 1, stored_matrix_element_count
                    call log_bits("H",0,symmetry_block_index,0,matrix_entry,hamiltonian_values(matrix_entry))
                end do
            end if
            call profile_stage("matrix_assembly",0,symmetry_block_index,started)
            ! Determine average energy
            hamiltonian_energy_shift = 0.0_real64
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(stored_segment_end_offsets(csf_index))
            end do
            started = mcp_now()
            call mpi_barrier(mpi_comm_world,mpi_error_code)
            call profile_stage("matrix_barrier",0,symmetry_block_index,started)
            started = mcp_now()
            call mpi_allreduce (hamiltonian_energy_shift,tmp,1,mpi_double_precision,mpi_sum,  &
                mpi_comm_world,mpi_error_code)
            call profile_stage("matrix_reduction",0,symmetry_block_index,started)
            hamiltonian_energy_shift = tmp / csf_count
            ! EAV = EAV/NCF
            mean_energy_by_block(symmetry_block_index) = hamiltonian_energy_shift

            !  Print Hamiltonian matrix and average energy
            !  hmout is not general

            if (mpi_rank == 0) write (*, '(" Average energy = ",1p,d18.10," Hartrees")') hamiltonian_energy_shift

            ! Subtract the average energy from the diagonal elements
            ! to reduce the condition number of the matrix
            !            EMT(IENDC(NPROCS*(I-1)+MYID+1)) = EMT(IENDC(NPROCS*(I-1)+MYID+1))&
            !                - EAV
            do i = mpi_rank + 1, csf_count, mpi_process_count
                idiag = stored_segment_end_offsets(i)     ! new mode: each row ends in diagonal
                hamiltonian_values(idiag) = hamiltonian_values(idiag) - hamiltonian_energy_shift
            end do

            ! Compute and store eigenpairs
            call diagonalize_hamiltonian_block_mpi (first_davidson_solve, ldbpg(3),                      &
                symmetry_block_index, block_csf_offset, block_eigenvalue_offset, block_mixing_offset, total_csf_count)

            ! Damp and Schmidt orthogonalise eigenvectors for OL calculations

            if (.not.first) then

                do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)

                    iofset = (level_index - 1)*csf_count + block_mixing_offset
                    jother = level_index

                    ! cdamp has the original non-block feature
                    cdampj = mixing_damping_factors(level_index + block_eigenvalue_offset)
                    if (cdampj == 0.0_real64) cycle         ! So SURE ???

                    omcdaj = 1.0_real64 - cdampj

                    ! ...Damp eigenvector and determine the new dominant component
                    do
                        amax = 0.0_real64
                        do i = 1, csf_count
                            evecij = omcdaj*csf_mixing_coefficients(i+iofset) + cdampj*cmvl(i + iofset)
                            csf_mixing_coefficients(i+iofset) = evecij
                            wa = abs(evecij)
                            if (wa <= amax) cycle
                            amax = wa
                            ia = i
                        end do

                        ! ...compute the normalization factor
                        sum = 0.0_real64
                        do i = 1, csf_count
                            sum = sum + csf_mixing_coefficients(i+iofset)**2
                        end do
                        dnfac = 1.0_real64/sqrt(sum)

                        ! ...Renormalize and invert as necessary
                        if (csf_mixing_coefficients(ia+iofset) < 0.0_real64) dnfac = -dnfac
                        call dscal (csf_count, dnfac, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)

                        ! ...Schmidt orthogonalise
                        do
                            jother = jother - 1
                            if (jother < 1) exit
                            jofset = (jother - 1)*csf_count + block_mixing_offset
                            ovrlap = ddot(csf_count,csf_mixing_coefficients(iofset+1:iofset+(csf_count)),1, &
                                & csf_mixing_coefficients(jofset+1:jofset+(csf_count)),1)
                            if (ovrlap /= 0.0_real64) then          ! So SURE ???
                                omcdaj = 1.0_real64
                                cdampj = -ovrlap
                                call dcopy (csf_count, csf_mixing_coefficients(jofset+1:jofset+(csf_count)), 1, cmvl(iofset + &
                                    & 1:iofset+(csf_count)), 1)
                                exit
                            end if
                        end do
                        if (jother < 1) exit
                    end do
                end do
            end if

            ! Write out the eigenpair information: ASF symmetries, eigenvalues,
            ! and eigenvectors to GRASP92 mixing coefficients File

            if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
                iattmp = 999
                iastmp = 999
            else
                iattmp = iabs(packed_block_symmetries(symmetry_block_index))
                if(packed_block_symmetries(symmetry_block_index) >= 0) then
                    iastmp = 1
                else
                    iastmp = -1
                end if
                ! IATTMP = IATJPO(NCMINPAT + 1)
                ! IASTMP = IASPAR(NCMINPAT + 1)
            end if

            if (mpi_rank == 0) then
                write (mixing_output_unit) symmetry_block_index, csf_count, eigenvector_counts_by_block(symmetry_block_index), &
                    & iattmp, iastmp
                write (mixing_output_unit) (requested_level_indices(i + block_eigenvalue_offset),i=1, &
                    eigenvector_counts_by_block(symmetry_block_index))
                write (mixing_output_unit) hamiltonian_energy_shift, (shifted_eigenvalues(i + block_eigenvalue_offset),i=1, &
                    eigenvector_counts_by_block(symmetry_block_index))
                write (mixing_output_unit) ((csf_mixing_coefficients(i+(level_index-1)*csf_count+block_mixing_offset),i=1, &
                    csf_count), &
                    & level_index=1, &
                    eigenvector_counts_by_block(symmetry_block_index)&
                    )
            end if

            call dalloc (hamiltonian_values, "hamiltonian_values", "solve_hamiltonian_blocks_mpi")
            call dalloc (stored_segment_end_offsets, "stored_segment_end_offsets", "solve_hamiltonian_blocks_mpi")
            call dalloc (matrix_partner_indices, "matrix_partner_indices", "solve_hamiltonian_blocks_mpi")

        end do
        ! Deallocate the temporary storage
        if (.not.first) then
            call dalloc (cmvl, "CMVL", "solve_hamiltonian_blocks_mpi")
        else
            first = .false.
        end if

        return
    end subroutine solve_hamiltonian_blocks_mpi
end module solve_hamiltonian_blocks_mpi_interface
