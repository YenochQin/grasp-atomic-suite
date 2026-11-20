module solve_hamiltonian_blocks_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_hamiltonian_blocks
contains

    !***********************************************************************
    !                                                                      *
    subroutine solve_hamiltonian_blocks(first_davidson_solve)
        ! Module dependencies
        use mcp_file_metadata, only: mcp_input_units
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
        !   Midified by G. Gaigalas                              05 Feb 2017   *
        !      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
        !                                  JCUPA(NNNW*NCF)                     *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use array_memory_management, only: alloc, dalloc
        use scf_damping_factors, only: mixing_damping_factors
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use debug_output_settings, only: ldbpg
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, symmetry_block_count
        use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
            matrix_partner_indices
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use serial_mpi_context, only: mpi_process_count, mpi_rank
        use orbital_csf_basis, only: csf_count, orbital_count
        use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, csf_offsets_by_block, &
            eigenvalue_offsets_by_block, mixing_offsets_by_block
        use block_mean_energies, only: mean_energy_by_block
        use level_symmetry_metadata,          only: packed_block_symmetries

        ! Procedure interfaces
        use assemble_hamiltonian_block_interface, only: assemble_hamiltonian_block
        use diagonalize_hamiltonian_block_interface, only: diagonalize_hamiltonian_block
        use blas_lapack_interfaces, only: dcopy, ddot, dscal
        implicit none (type, external)

        ! Arguments
        logical, intent(in) ::  first_davidson_solve

        ! Local declarations
        integer :: mcp_channel, file_electron_count, discarded_header_value, symmetry_block_index, file_block_index, i, &
            & block_csf_offset, block_eigenvalue_offset, &
            block_mixing_offset, file_csf_count, file_coefficient_count, lab, contribution_count, itmp, csf_index, idiag, &
                & level_index, iofset, &
            jother, ia, jofset, iattmp, iastmp, nelmntgg
        real(real64), pointer :: cmvl(:)
        real(real64) :: amax, cdampj, dnfac, evecij, sum,  tmp, omcdaj, ovrlap, wa
        logical :: first
        character :: mcp_channel_label*3

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
            call alloc (cmvl, mixing_coefficient_storage_size, "CMVL", "MATRIX")
            call dcopy (mixing_coefficient_storage_size, csf_mixing_coefficients, 1, cmvl, 1)
        end if

        !=======================================================================
        !   Position the files - MCP files (unit NFILE) for reading
        !   and mixing coefficients file (unit 25) for writing
        !=======================================================================

        do mcp_channel = 30, 32 + max_mcp_file_multipole_rank
            rewind (mcp_input_units(mcp_channel))
            if (mcp_channel == 30) then
                read (mcp_input_units(mcp_channel))
                read (mcp_input_units(mcp_channel))
                read (mcp_input_units(mcp_channel))
            end if
            read (mcp_input_units(mcp_channel))
            read (mcp_input_units(mcp_channel))
            read (mcp_input_units(mcp_channel))
        end do

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

        !------------------------------------------------
        do symmetry_block_index = 1, symmetry_block_count                      ! block do-loop
            !------------------------------------------------

            !=======================================================================
            !   Read indeces of non-zero elements from mcp.30 file. Note the
            !   format has been changed to lower-triangle-by-rows.
            ! Length of iendc can be reduced
            !=======================================================================

            read (mcp_input_units(30)) mcp_channel_label, file_block_index, csf_count
            if (file_block_index/=symmetry_block_index .or. csf_count/=csf_counts_by_block(symmetry_block_index)) error stop &
                "matrx: jblockt .NE. jblock .OR. ncf1 .NE. ncf2"
            read (mcp_input_units(30)) nelmntgg
            stored_matrix_element_count = int(nelmntgg,kind=int64)
            call alloc (matrix_partner_indices, stored_matrix_element_count, "IROW", "MATRIX")
            call alloc (hamiltonian_values, stored_matrix_element_count, "EMT", "MATRIX")
            !cjb          ALLOC (IENDC, NCF + 1,...)  ->  ALLOC (IENDC, 0, NCF,...)
            !cjb     CALL ALLOC (IENDC, NCF + 1, 'IENDC', 'MATRIX' )
            call alloc (stored_segment_end_offsets, 0, csf_count, "IENDC", "MATRIX" )

            !        ! may not be necessary if iendc is ALWAYS used
            stored_segment_end_offsets(0:csf_count) = 0                        ! the way it is assigned here.

            !...EMT will be accumulated in setham
            hamiltonian_values(:stored_matrix_element_count) = 0.0_real64

            read (mcp_input_units(30)) (stored_segment_end_offsets(i),i=mpi_rank + 1,csf_count,mpi_process_count), &
                & (matrix_partner_indices(i),i=1,stored_matrix_element_count)

            block_csf_offset = csf_offsets_by_block(symmetry_block_index)
            block_eigenvalue_offset = eigenvalue_offsets_by_block(symmetry_block_index)
            block_mixing_offset = mixing_offsets_by_block(symmetry_block_index)

            !=======================================================================
            !   Skip current block if no eigenlaue is required
            !=======================================================================

            if (eigenvector_counts_by_block(symmetry_block_index) == 0) then
                do mcp_channel = 31, 32 + max_mcp_file_multipole_rank
                    read (mcp_input_units(mcp_channel)) mcp_channel_label, file_block_index, file_csf_count, file_coefficient_count
                    if (file_block_index /= symmetry_block_index) error stop "matrx: jblockt .NE. jblock"
                    if (file_csf_count /= csf_count) error stop "matrx: ncft .NE. ncf"

                    read (mcp_input_units(mcp_channel)) lab, contribution_count
                    do while(lab/=0 .or. contribution_count/=0)
                        read (mcp_input_units(mcp_channel)) (itmp,itmp,tmp,i=1,contribution_count)
                        read (mcp_input_units(mcp_channel)) lab, contribution_count
                    end do
                end do

                call dalloc (stored_segment_end_offsets, "IENDC", "MATRIX")
                call dalloc (matrix_partner_indices, "IROW", "MATRIX")

                cycle

            end if

            !=======================================================================
            !   Generate the Hamiltonian matrix - average energy is removed here
            !=======================================================================

            call assemble_hamiltonian_block (symmetry_block_index, mpi_rank, mpi_process_count)
            !
            !   Determine average energy
            !
            hamiltonian_energy_shift = 0.0_real64
            do csf_index = mpi_rank + 1, csf_count, mpi_process_count
                hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(stored_segment_end_offsets(csf_index))
            end do
            !         DO IR = 1, (NCF - (MYID + 1) + NPROCS)/NPROCS
            !            EAV = EAV + EMT(IENDC(NPROCS*(IR-1)+MYID+1))
            !         END DO

            hamiltonian_energy_shift = hamiltonian_energy_shift/csf_count
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

            call diagonalize_hamiltonian_block (first_davidson_solve, ldbpg(3),                      &
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
                    !cjb iastmp error corrected by Pawel Syty
                    !PS            iastmp = 2
                    iastmp = 1
                else
                    !cjb iastmp error corrected by Pawel Syty
                    !PS            iastmp = 1
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

            call dalloc (hamiltonian_values, "EMT", "MATRIX")
            call dalloc (stored_segment_end_offsets, "IENDC", "MATRIX")
            call dalloc (matrix_partner_indices, "IROW", "MATRIX")

            !----------------------
        end do
        !----------------------
        !
        !   Deallocate the temporary storage
        !
        if (.not.first) then
            call dalloc (cmvl, "CMVL", "MATRIX")
        else
            first = .false.
        end if


        return
    end subroutine solve_hamiltonian_blocks
end module solve_hamiltonian_blocks_interface
