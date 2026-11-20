! Stores dimensions and selected eigenvector counts for the Hamiltonian symmetry blocks.
module hamiltonian_symmetry_blocks
    implicit none (type, external)

    ! Module state
    integer :: symmetry_block_count
    integer, pointer :: csf_counts_by_block(:) => null()
    integer, pointer :: eigenvector_counts_by_block(:), highest_selected_level_by_block(:)
end module hamiltonian_symmetry_blocks
