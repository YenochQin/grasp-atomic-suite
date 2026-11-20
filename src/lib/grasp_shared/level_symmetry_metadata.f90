! Stores 2J+1 and parity signs for the active, final, and initial sets of levels.
module level_symmetry_metadata
    implicit none (type, external)

    ! Module constants
    integer, parameter::max_symmetry_blocks = 50  ! Capacity of the symmetry-block lookup array.

    ! Module state
    integer, pointer :: level_two_j_plus_one(:), level_parity_signs(:)
    integer, pointer :: final_level_two_j_plus_one(:), final_level_parity_signs(:)
    integer, pointer :: initial_level_two_j_plus_one(:), initial_level_parity_signs(:)
    integer :: packed_block_symmetries(max_symmetry_blocks)
end module level_symmetry_metadata
