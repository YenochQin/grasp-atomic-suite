module solve_hamiltonian_blocks_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_hamiltonian_blocks_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine solve_hamiltonian_blocks_mpi(first_davidson_solve)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use calculation_file_units, only: mixing_output_unit
        !                                                                      *
        !   Calls routines to form the Hamiltonian matrix and to diagonalise   *
        !   it. The total angular momenta and parity of each  ASF  is found;   *
        !   the eigenvectors are normalised so  that the sign of the largest   *
        !   element is positive.                                               *
        !                                                                      *
        !   Call(s) to: [LIB92] ALLOC, DALLOC.                                 *
        !               [RSCF92]: SETHAM, HMOUT.                               *
        !               [BLAS]: DCOPY/SCOPY, DSCAL/SSCAL                       *
        !                                                                      *
        !                                         Last revision: 24 Dec 1992   *
        ! Block version by Xinghong He                           07 Aug 1998   *
        !   Modified by G. Gaigalas               Last revision: 07 Sep 2021   *
        !      Spin-angular coefficiants have been put on memory               *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use array_memory_management, only: alloc, dalloc
        use scf_damping_factors, only: mixing_damping_factors
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use debug_output_settings, only: ldbpg
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use hamiltonian_symmetry_blocks, only: eigenvector_counts_by_block, symmetry_block_count
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
            matrix_partner_indices
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use mpi_state, only: mpi_allreduce, mpi_barrier, mpi_comm_world, mpi_double_precision, mpi_error_code, &
            mpi_process_count, mpi_rank, mpi_sum
        use orbital_csf_basis, only: csf_count, orbital_count
        use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, csf_offsets_by_block, &
            eigenvalue_offsets_by_block, mixing_offsets_by_block
        use block_mean_energies, only: mean_energy_by_block
        use level_symmetry_metadata,          only: packed_block_symmetries
        use rmcdhf_mcp_memory,    only: imin1_30, imax1_30,imin2_30, imax2_30

        ! Procedure interfaces
        use assemble_hamiltonian_block_interface, only: assemble_hamiltonian_block
        use diagonalize_hamiltonian_block_mpi_interface, only: diagonalize_hamiltonian_block_mpi
        use load_mcp_records_into_memory_interface, only: load_mcp_records_into_memory
        use read_mcp_coefficients_from_memory_interface, only: read_mcp_coefficients_from_memory
        use blas_lapack_interfaces, only: dcopy, ddot, dscal
        implicit none (type, external)

        ! Arguments
        logical, intent(in) ::  first_davidson_solve

        ! Local declarations
        integer :: mcp_channel, file_electron_count, discarded_header_value, symmetry_block_index, i, block_csf_offset, &
            block_eigenvalue_offset, block_mixing_offset, lab, contribution_count, csf_index, idiag, level_index, iofset, &
            jother, ia, jofset, iattmp, iastmp, nelmntgg
        real(real64), pointer :: cmvl(:)
        real(real64) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
        logical :: first
        integer(int64) :: total_contribution_count(20), total_secondary_contribution_count(20)

        save first

        ! Initial values
        !
        data first/ .true./
        !
        !     POINTER (cmvl(1))
        !
        !-----------------------------------------------------------------------

        if (mpi_rank == 0) write (output_unit, *)

        !   Allocate memory for CMVL once (the maximum size)
        !   Save previous estimate of eigenvectors

        if (.not.first) then
            call alloc (cmvl, mixing_coefficient_storage_size, "CMVL", "solve_hamiltonian_blocks_mpi")
            call dcopy (mixing_coefficient_storage_size, csf_mixing_coefficients, 1, cmvl, 1)
        end if

        ! To put in ncmin and nvecsiz. Values read here are the same
        ! as those from elsewhere (shch as common blocks)
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

        !=======================================================================
        !   Do the job block by block
        !=======================================================================

        total_contribution_count = 0
        total_secondary_contribution_count = 0
        !------------------------------------------------
        do symmetry_block_index = 1, symmetry_block_count                      ! block do-loop
            !------------------------------------------------

            !=======================================================================
            !   Read indeces of non-zero elements from mcp.30 file. Note the
            !   format has been changed to lower-triangle-by-rows.
            ! Length of iendc can be reduced
            !=======================================================================

            nelmntgg = imax2_30(symmetry_block_index) - imin2_30(symmetry_block_index)
            stored_matrix_element_count = int(nelmntgg,kind=int64)
            csf_count = imax1_30(symmetry_block_index) - imin1_30(symmetry_block_index)
            call alloc (matrix_partner_indices, stored_matrix_element_count, "matrix_partner_indices", &
                "solve_hamiltonian_blocks_mpi")
            call alloc (hamiltonian_values, stored_matrix_element_count, "hamiltonian_values", "solve_hamiltonian_blocks_mpi")
            call alloc (stored_segment_end_offsets, 0, csf_count, "stored_segment_end_offsets", "solve_hamiltonian_blocks_mpi" )

            !                                                ! may not be necessary if iendc is ALWAYS used
            !GG         IENDC(:NCF) = 0                        ! the way it is assigned here.
            stored_segment_end_offsets(0:csf_count) = 0

            !...EMT will be accumulated in setham
            hamiltonian_values(:stored_matrix_element_count) = 0.0_real64
            call read_mcp_coefficients_from_memory(symmetry_block_index)

            block_csf_offset = csf_offsets_by_block(symmetry_block_index)
            block_eigenvalue_offset = eigenvalue_offsets_by_block(symmetry_block_index)
            block_mixing_offset = mixing_offsets_by_block(symmetry_block_index)

            !=======================================================================
            !   Skip current block if no eigenlaue is required
            !=======================================================================

            if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
                do mcp_channel = 31, 32 + max_mcp_file_multipole_rank
                    total_secondary_contribution_count(mcp_channel-30) = total_secondary_contribution_count(mcp_channel-30) + 1
                    call load_mcp_records_into_memory(mcp_channel, total_secondary_contribution_count(mcp_channel-30), lab, &
                        & contribution_count)
                    do while(lab/=0 .or. contribution_count/=0)
                        total_contribution_count(mcp_channel-30) = total_contribution_count(mcp_channel-30) + contribution_count
                        total_secondary_contribution_count(mcp_channel-30) = total_secondary_contribution_count(mcp_channel-30) + 1
                        call load_mcp_records_into_memory(mcp_channel, total_secondary_contribution_count(mcp_channel-30), &
                            & lab, contribution_count)
                    end do
                end do

                call dalloc (stored_segment_end_offsets, "stored_segment_end_offsets", "solve_hamiltonian_blocks_mpi")
                call dalloc (matrix_partner_indices, "matrix_partner_indices", "solve_hamiltonian_blocks_mpi")

                cycle

            end if

            !=======================================================================
            !   Generate the Hamiltonian matrix - average energy is removed here
            !=======================================================================

            call assemble_hamiltonian_block (symmetry_block_index, mpi_rank, mpi_process_count, total_contribution_count, &
                & total_secondary_contribution_count)
            !
            !   Determine average energy
            !
            hamiltonian_energy_shift = 0.0_real64
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(stored_segment_end_offsets(csf_index))
            end do
            call mpi_barrier(mpi_comm_world,mpi_error_code)
            call mpi_allreduce (hamiltonian_energy_shift,tmp,1,mpi_double_precision,mpi_sum,  &
                mpi_comm_world,mpi_error_code)
            hamiltonian_energy_shift = tmp / csf_count
            !GG         EAV = EAV/NCF
            mean_energy_by_block(symmetry_block_index) = hamiltonian_energy_shift

            ! Print Hamiltonian matrix and average energy
            ! hmout is not general
            !call hmout (0, 1, ncf)

            if (mpi_rank == 0) write (*, '(" Average energy = ",1p,d18.10," Hartrees")') hamiltonian_energy_shift

            ! Subtract the average energy from the diagonal elements
            ! to reduce the condition number of the matrix
            !         DO I = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
            !            EMT(IENDC(NPROCS*(I-1)+MYID+1)) = EMT(IENDC(NPROCS*(I-1)+MYID+1))&
            !                - EAV
            do i = mpi_rank + 1, csf_count, mpi_process_count
                idiag = stored_segment_end_offsets(i)     ! new mode: each row ends in diagonal
                hamiltonian_values(idiag) = hamiltonian_values(idiag) - hamiltonian_energy_shift
            end do

            !=======================================================================
            !   Compute and store eigenpairs
            !=======================================================================
            call diagonalize_hamiltonian_block_mpi (first_davidson_solve, ldbpg(3),                      &
                symmetry_block_index, block_csf_offset, block_eigenvalue_offset, block_mixing_offset, total_csf_count)

            !=======================================================================
            !   Damp and Schmidt orthogonalise eigenvectors for OL calculations
            !=======================================================================

            if (.not.first) then

                do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)

                    iofset = (level_index - 1)*csf_count + block_mixing_offset
                    jother = level_index

                    ! cdamp has the original non-block feature
                    cdampj = mixing_damping_factors(level_index + block_eigenvalue_offset)
                    if (cdampj == 0.0_real64) cycle         ! So SURE ???

                    omcdaj = 1.0_real64 - cdampj

                    !...Damp eigenvector and determine the new dominant component
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

                        !...compute the normalization factor
                        sum = 0.0_real64
                        do i = 1, csf_count
                            sum = sum + csf_mixing_coefficients(i+iofset)**2
                        end do
                        dnfac = 1.0_real64/sqrt(sum)

                        !...Renormalize and invert as necessary
                        if (csf_mixing_coefficients(ia+iofset) < 0.0_real64) dnfac = -dnfac
                        call dscal (csf_count, dnfac, csf_mixing_coefficients(iofset+1:iofset+(csf_count)), 1)

                        !...Schmidt orthogonalise
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

            !   Write out the eigenpair information: ASF symmetries, eigenvalues,
            !   and eigenvectors to GRASP92 mixing coefficients File

            if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
                iattmp = 999
                iastmp = 999
            else
                !GGGG
                iattmp = iabs(packed_block_symmetries(symmetry_block_index))
                if(packed_block_symmetries(symmetry_block_index) >= 0) then
                    iastmp = 1
                else
                    iastmp = -1
                end if
                !GG            IATTMP = IATJPO(NCMINPAT + 1)
                !GG            IASTMP = IASPAR(NCMINPAT + 1)
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

            !----------------------
        end do
        !----------------------
        !
        !   Deallocate the temporary storage
        !
        if (.not.first) then
            call dalloc (cmvl, "CMVL", "solve_hamiltonian_blocks_mpi")
        else
            first = .false.
        end if


        return
    end subroutine solve_hamiltonian_blocks_mpi
end module solve_hamiltonian_blocks_mpi_interface
