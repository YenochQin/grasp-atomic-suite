module resize_transition_coefficient_arrays_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: resize_transition_coefficient_arrays
contains

    ! This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),
    ! and deallocates (IMODE = 3) storage for  certain arrays that are
    ! local to SUBROUTINE TRSORT.
    !
    ! Written by Farid A. Parpia. Last revision: 28 Dec 1992
    subroutine resize_transition_coefficient_arrays(coefficient_sort_links, initial_csf_index_workspace, final_csf_indices, &
        & sorted_angular_coefficients, coefficient_capacity, resize_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use array_memory_management, only: alloc, dalloc, ralloc
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: coefficient_capacity
        integer , intent(in) :: resize_mode
        integer, pointer :: coefficient_sort_links(:), initial_csf_index_workspace(:), final_csf_indices(:)
        real(real64), pointer :: sorted_angular_coefficients(:)

        ! Local declarations
        integer :: new_capacity

        select case (resize_mode)
        case (1)
            ! Initial array dimension
            coefficient_capacity = 64
            ! Allocate storage for arrays
            call alloc (coefficient_sort_links, coefficient_capacity, "coefficient_sort_links", &
                "resize_transition_coefficient_arrays")
            call alloc (initial_csf_index_workspace, coefficient_capacity, "initial_csf_index_workspace", &
                "resize_transition_coefficient_arrays")
            call alloc (final_csf_indices, coefficient_capacity, "final_csf_indices", "resize_transition_coefficient_arrays")
            call alloc (sorted_angular_coefficients, coefficient_capacity, "sorted_angular_coefficients", &
                "resize_transition_coefficient_arrays")
        case (2)
            ! Double the allocation of storage for the arrays
            new_capacity = 2*coefficient_capacity
            call ralloc (coefficient_sort_links, new_capacity, "coefficient_sort_links", "resize_transition_coefficient_arrays")
            call ralloc (initial_csf_index_workspace, new_capacity, "initial_csf_index_workspace", &
                "resize_transition_coefficient_arrays")
            call ralloc (final_csf_indices, new_capacity, "final_csf_indices", "resize_transition_coefficient_arrays")
            call ralloc (sorted_angular_coefficients, new_capacity, "sorted_angular_coefficients", &
                "resize_transition_coefficient_arrays")
            coefficient_capacity = new_capacity
        case (3)
            ! Deallocate the storage for the arrays
            call dalloc (coefficient_sort_links,  "coefficient_sort_links", "resize_transition_coefficient_arrays")
            call dalloc (initial_csf_index_workspace, "initial_csf_index_workspace", "resize_transition_coefficient_arrays")
            call dalloc (final_csf_indices, "final_csf_indices", "resize_transition_coefficient_arrays")
            call dalloc (sorted_angular_coefficients,  "sorted_angular_coefficients", "resize_transition_coefficient_arrays")
        case default
            write (output_unit, *) "resize_transition_coefficient_arrays: Invalid resize_mode = ", resize_mode
            stop
        end select
        return
    end subroutine resize_transition_coefficient_arrays
end module resize_transition_coefficient_arrays_interface
