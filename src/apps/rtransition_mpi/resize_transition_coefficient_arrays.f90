module resize_transition_coefficient_arrays_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: resize_transition_coefficient_arrays
contains

    !***********************************************************************
    !                                                                      *
    subroutine resize_transition_coefficient_arrays(coefficient_sort_links, initial_csf_index_workspace, final_csf_indices, &
        & sorted_angular_coefficients, coefficient_capacity, resize_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
        !   and deallocates (IMODE = 3) storage for  certain arrays that are   *
        !   local to SUBROUTINE TRSORT.                                        *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
        !                                                                      *
        !***********************************************************************
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
        !-----------------------------------------------
        !
        select case (resize_mode)
        case (1)
            !
            !   Initial array dimension
            !
            coefficient_capacity = 64
            !
            !   Allocate storage for arrays
            !
            call alloc (coefficient_sort_links, coefficient_capacity, "IPTR", "ALCNMA")
            call alloc (initial_csf_index_workspace, coefficient_capacity, "ISLDR", "ALCNMA")
            call alloc (final_csf_indices, coefficient_capacity, "ISLDR1", "ALCNMA")
            call alloc (sorted_angular_coefficients, coefficient_capacity, "XSLDR", "ALCNMA")
        !
        case (2)
            !
            !   Double the allocation of storage for the arrays
            !
            new_capacity = 2*coefficient_capacity
            !
            call ralloc (coefficient_sort_links, new_capacity, "IPTR", "ALCNMA")
            call ralloc (initial_csf_index_workspace, new_capacity, "ISLDR", "ALCNMA")
            call ralloc (final_csf_indices, new_capacity, "ISLDR1", "ALCNMA")
            call ralloc (sorted_angular_coefficients, new_capacity, "XSLDR", "ALCNMA")
            !
            coefficient_capacity = new_capacity
        !
        case (3)
            !
            !   Deallocate the storage for the arrays
            !
            call dalloc (coefficient_sort_links,  "IPTR", "ALCNMA")
            call dalloc (initial_csf_index_workspace, "ISLDR", "ALCNMA")
            call dalloc (final_csf_indices, "ISLDR1", "ALCNMA")
            call dalloc (sorted_angular_coefficients,  "XSLDR", "ALCNMA")
        !
        case default
            !
            write (output_unit, *) "ALCNMA: Invalid argument IMODE = ", resize_mode
            stop
            !
        end select
        !
        return
    end subroutine resize_transition_coefficient_arrays
end module resize_transition_coefficient_arrays_interface
