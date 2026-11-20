module read_ci_mixing_coefficients_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_mixing_coefficients_mpi
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_ci_mixing_coefficients_mpi (block_labels)
        ! Module dependencies
        use calculation_file_units, only: mixing_output_unit
        !                                                                      *
        !   Determines the eigenpairs required;  this information is written   *
        !   to the head of the  .mix  file.                                    *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, lodstate                               *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 24 Nov 1992   *
        !   Block version by Xinghong He          Last revision:  9 Jun 1998   *
        !                                                                      *
        !***********************************************************************
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: electron_count, requested_level_indices, selected_level_count
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use transition_basis_blocks,           only: level_block_indices
        use orbital_csf_basis, only: orbital_count
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_integer4, mpi_rank
        use block_level_selection, only: read_block_level_selection
        implicit none (type, external)

        ! Arguments
        character(len=8), intent(in):: block_labels(*)

        ! Local declarations
        integer :: total_csf_count, noffset, vector_size, symmetry_block_index, level_index
        integer :: nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
        !-----------------------------------------------

        ! lodstate generates
        !    nevblk(), ncmaxblk()
        !    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
        ! Thus we let node-0 do it and then broadcast here

        call alloc (highest_selected_level_by_block, symmetry_block_count, "NCMAXBLK", "LODMIXmpi")
        call alloc (eigenvector_counts_by_block, symmetry_block_count, "NEVBLK", "LODMIXmpi")

        !      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
        !      stop
        if (mpi_rank == 0) then
            call read_block_level_selection (block_labels)
        end if

        call mpi_bcast (eigenvector_counts_by_block, symmetry_block_count, mpi_integer, 0,               &
            mpi_comm_world, mpi_error_code)
        call mpi_bcast (highest_selected_level_by_block, symmetry_block_count, mpi_integer, 0,             &
            mpi_comm_world, mpi_error_code)
        call mpi_bcast (selected_level_count, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)

        if (selected_level_count /= 0) then
            if (mpi_rank /= 0) then   ! lodstate allocated it from node-0
                call alloc (requested_level_indices, selected_level_count, "ICCMIN", "LODMIXmpi")
            end if

            call mpi_bcast (requested_level_indices, selected_level_count, mpi_integer4, 0,            &
                mpi_comm_world, mpi_error_code)
        end if

        !
        !   Determine other auxiliary arrays, parameters
        !
        !    idxblk() is the block number of an eigenstate
        !    nvecsiz  is the total size of the eigenvector array
        !
        call alloc (level_block_indices, selected_level_count, "IDXBLK", "LODMIXmpi")
        total_csf_count = 0
        noffset = 0
        vector_size = 0
        do symmetry_block_index = 1, symmetry_block_count
            do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
                level_block_indices(level_index + noffset) = symmetry_block_index
            end do
            total_csf_count = total_csf_count + csf_counts_by_block(symmetry_block_index)
            noffset = noffset + eigenvector_counts_by_block(symmetry_block_index)
            vector_size = vector_size + csf_counts_by_block(symmetry_block_index) * eigenvector_counts_by_block( &
                & symmetry_block_index)
        end do
        if (noffset /= selected_level_count) stop "lodmixmpi: ncmin trouble"
        !
        !   Write header data to the  .mix  file
        !
        if (mpi_rank == 0) then
            nelec4 = electron_count
            ncftot4 = total_csf_count
            nw4 = orbital_count
            ncmin4 = selected_level_count
            nvecsiz4 = vector_size
            nblock4 = symmetry_block_count
            write (mixing_output_unit) nelec4, ncftot4, nw4, ncmin4, nvecsiz4, nblock4
        end if
        !
        return
    end subroutine read_ci_mixing_coefficients_mpi
end module read_ci_mixing_coefficients_mpi_interface
