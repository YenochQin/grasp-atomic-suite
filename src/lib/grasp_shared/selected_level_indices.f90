! Stores the selected level counts and indices for the active, final, and initial bases.
module selected_level_indices
    implicit none (type, external)

    ! Module state
    integer :: active_level_count, highest_active_level_index
    integer :: final_level_count, highest_final_level_index
    integer :: initial_level_count, highest_initial_level_index
    integer, pointer :: active_level_indices(:), final_level_indices(:), initial_level_indices(:)
end module selected_level_indices
