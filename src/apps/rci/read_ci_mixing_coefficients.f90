module read_ci_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_mixing_coefficients
contains

    !***********************************************************************
    !                                                                      *
    subroutine read_ci_mixing_coefficients (block_labels)
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
        use atomic_calculation_state, only: electron_count, selected_level_count
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use transition_basis_blocks, only: level_block_indices
        use orbital_csf_basis, only: orbital_count
        use block_level_selection, only: read_block_level_selection
        implicit none (type, external)

        ! Arguments
        character(len=8), intent(in):: block_labels(*)

        ! Local declarations
        integer :: total_csf_count, noffset, vector_size, symmetry_block_index, level_index
        !-----------------------------------------------

        ! lodstate generates
        !    nevblk(), ncmaxblk()
        !    ncmin, iccmin(1:ncmin) -- via items (memories allocated there)
        ! Thus we let node-0 do it and then broadcast here

        call alloc (highest_selected_level_by_block, symmetry_block_count, "highest_selected_level_by_block", &
            "read_ci_mixing_coefficients")
        call alloc (eigenvector_counts_by_block, symmetry_block_count, "eigenvector_counts_by_block", "read_ci_mixing_coefficients")

        !      print *, ' LODMIX: change lodstate arg-list - see rscf2/getold.f90 '
        !      stop
        call read_block_level_selection (block_labels)
        !     CALL lodstate (nblock, ncfblk(1), idblk, nevblk, ncmaxblk)
        !
        !   Determine other auxiliary arrays, parameters
        !
        !    idxblk() is the block number of an eigenstate
        !    nvecsiz  is the total size of the eigenvector array
        !
        call alloc (level_block_indices, selected_level_count, "level_block_indices", "read_ci_mixing_coefficients")
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
        if (noffset /= selected_level_count) stop "read_ci_mixing_coefficients: selected level count mismatch"
        !
        !   Write header data to the  .mix  file
        !
        write (mixing_output_unit) electron_count, total_csf_count, orbital_count, selected_level_count, vector_size, &
            symmetry_block_count
        !
        return
    end subroutine read_ci_mixing_coefficients
end module read_ci_mixing_coefficients_interface
