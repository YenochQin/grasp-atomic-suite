! Stores CSF and eigenvector counts, symmetries, and index arrays for transition bases.
module transition_basis_blocks
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: cumulative_csf_counts(20), transition_level_counts_by_block(20), transition_csf_counts_by_block(20), &
        transition_two_j_by_block(20)
    integer :: transition_block_count
    real(real64), pointer :: level_block_indices(:)
    integer :: mixing_file_electron_count, mixing_file_csf_count, mixing_file_orbital_count, mixing_file_level_count, &
        mixing_file_block_count
    character(64) :: transition_csf_basis_ids(2) = ""
    integer :: initial_transition_block_count
    integer :: final_transition_block_count
    integer :: initial_cumulative_csf_counts(10)
    integer :: final_cumulative_csf_counts(10)
end module transition_basis_blocks
