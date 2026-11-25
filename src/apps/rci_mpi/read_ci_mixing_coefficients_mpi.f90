module read_ci_mixing_coefficients_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_mixing_coefficients_mpi
contains

    ! Determines the eigenpairs required;  this information is written
    ! to the head of the  .mix  file.
    !
    ! Written by Farid A. Parpia. Last revision: 24 Nov 1992
    ! Block version by Xinghong He. Last revision:  9 Jun 1998
    subroutine read_ci_mixing_coefficients_mpi (block_labels, mixing_header)
        ! Module dependencies
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: electron_count, requested_level_indices, selected_level_count
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use transition_basis_blocks,           only: level_block_indices
        use orbital_csf_basis, only: orbital_count
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_integer4, mpi_rank
        use block_level_selection, only: read_block_level_selection
        use mpi_io_errors, only: abort_input_mpi
        implicit none (type, external)

        ! Arguments
        character(len=8), intent(in) :: block_labels(*)
        integer, intent(out) :: mixing_header(6)

        ! Local declarations
        integer :: total_csf_count, noffset, vector_size, symmetry_block_index, level_index

        ! lodstate generates
        !    nevblk(), ncmaxblk()
        !    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
        ! Thus we let node-0 do it and then broadcast here

        call alloc (highest_selected_level_by_block, symmetry_block_count, "highest_selected_level_by_block", &
            "read_ci_mixing_coefficients_mpi")
        call alloc (eigenvector_counts_by_block, symmetry_block_count, "eigenvector_counts_by_block", &
            "read_ci_mixing_coefficients_mpi")

        if (mpi_rank == 0) then
            call read_block_level_selection (block_labels, abort_input_mpi)
        end if

        call mpi_bcast (eigenvector_counts_by_block, symmetry_block_count, mpi_integer, 0,               &
            mpi_comm_world, mpi_error_code)
        call mpi_bcast (highest_selected_level_by_block, symmetry_block_count, mpi_integer, 0,             &
            mpi_comm_world, mpi_error_code)
        call mpi_bcast (selected_level_count, 1, mpi_integer, 0, mpi_comm_world, mpi_error_code)

        if (selected_level_count /= 0) then
            if (mpi_rank /= 0) then   ! lodstate allocated it from node-0
                call alloc (requested_level_indices, selected_level_count, "requested_level_indices", &
                    "read_ci_mixing_coefficients_mpi")
            end if

            call mpi_bcast (requested_level_indices, selected_level_count, mpi_integer4, 0,            &
                mpi_comm_world, mpi_error_code)
        end if

        ! Determine other auxiliary arrays, parameters
        !
        !  idxblk() is the block number of an eigenstate
        !  nvecsiz  is the total size of the eigenvector array
        call alloc (level_block_indices, selected_level_count, "level_block_indices", "read_ci_mixing_coefficients_mpi")
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
        if (noffset /= selected_level_count) stop "read_ci_mixing_coefficients_mpi: selected level count mismatch"
        mixing_header = [electron_count, total_csf_count, orbital_count, selected_level_count, vector_size, &
            symmetry_block_count]
        return
    end subroutine read_ci_mixing_coefficients_mpi
end module read_ci_mixing_coefficients_mpi_interface
