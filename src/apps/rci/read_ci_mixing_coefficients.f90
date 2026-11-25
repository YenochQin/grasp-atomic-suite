module read_ci_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_mixing_coefficients
contains

    ! Determine which eigenpairs RCI will calculate.
    !
    ! Written by Farid A. Parpia. Last revision: 24 Nov 1992
    ! Block version by Xinghong He. Last revision:  9 Jun 1998
    subroutine read_ci_mixing_coefficients (block_labels)
        ! Module dependencies
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: selected_level_count
        use hamiltonian_symmetry_blocks, only: eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use transition_basis_blocks, only: level_block_indices
        use block_level_selection, only: read_block_level_selection
        implicit none (type, external)

        ! Arguments
        character(len=8), intent(in) :: block_labels(*)

        ! Local declarations
        integer :: noffset, symmetry_block_index, level_index

        ! lodstate generates
        !    nevblk(), ncmaxblk()
        !    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
        ! Thus we let node-0 do it and then broadcast here

        call alloc (highest_selected_level_by_block, symmetry_block_count, "highest_selected_level_by_block", &
            "read_ci_mixing_coefficients")
        call alloc (eigenvector_counts_by_block, symmetry_block_count, "eigenvector_counts_by_block", "read_ci_mixing_coefficients")

        call read_block_level_selection (block_labels)
        ! Determine other auxiliary arrays, parameters
        !
        !  idxblk() is the block number of an eigenstate
        call alloc (level_block_indices, selected_level_count, "level_block_indices", "read_ci_mixing_coefficients")
        noffset = 0
        do symmetry_block_index = 1, symmetry_block_count
            do level_index = 1, eigenvector_counts_by_block(symmetry_block_index)
                level_block_indices(level_index + noffset) = symmetry_block_index
            end do
            noffset = noffset + eigenvector_counts_by_block(symmetry_block_index)
        end do
        if (noffset /= selected_level_count) stop "read_ci_mixing_coefficients: selected level count mismatch"
        return
    end subroutine read_ci_mixing_coefficients
end module read_ci_mixing_coefficients_interface
