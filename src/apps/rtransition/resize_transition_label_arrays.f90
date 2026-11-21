module resize_transition_label_arrays_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: resize_transition_label_arrays
contains

    !***********************************************************************
    !                                                                      *
    subroutine resize_transition_label_arrays(first_coefficient_indices, node_integral_labels, parent_nodes, &
        & left_child_nodes, coefficient_end_offsets, right_child_nodes, sorted_integral_labels, &
        label_capacity, resize_mode)
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
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: label_capacity
        integer , intent(in) :: resize_mode
        integer, pointer :: first_coefficient_indices(:), node_integral_labels(:), parent_nodes(:), left_child_nodes(:), &
            coefficient_end_offsets(:), right_child_nodes(:), sorted_integral_labels(:)

        ! Local declarations
        integer :: new_capacity
        !-----------------------------------------------
        !
        !
        select case (resize_mode)
        case (1)
            !
            !   Initial array dimension
            !
            label_capacity = 64
            !
            !   Allocate storage for arrays
            !
            call alloc (first_coefficient_indices, label_capacity, "first_coefficient_indices", "resize_transition_label_arrays")
            call alloc (node_integral_labels, label_capacity, "node_integral_labels", "resize_transition_label_arrays")
            call alloc (parent_nodes, label_capacity, "parent_nodes", "resize_transition_label_arrays")
            call alloc (left_child_nodes, label_capacity, "left_child_nodes", "resize_transition_label_arrays")
            call alloc (coefficient_end_offsets, label_capacity, "coefficient_end_offsets", "resize_transition_label_arrays")
            call alloc (right_child_nodes, label_capacity, "right_child_nodes", "resize_transition_label_arrays")
            call alloc (sorted_integral_labels, label_capacity, "sorted_integral_labels", "resize_transition_label_arrays")
        !
        case (2)
            !
            !   Double the allocation of storage for the arrays
            !
            new_capacity = 2*label_capacity
            !
            call ralloc (first_coefficient_indices, new_capacity, "first_coefficient_indices", "resize_transition_label_arrays")
            call ralloc (node_integral_labels, new_capacity, "node_integral_labels", "resize_transition_label_arrays")
            call ralloc (parent_nodes, new_capacity, "parent_nodes", "resize_transition_label_arrays")
            call ralloc (left_child_nodes, new_capacity, "left_child_nodes", "resize_transition_label_arrays")
            call ralloc (coefficient_end_offsets, new_capacity, "coefficient_end_offsets", "resize_transition_label_arrays")
            call ralloc (right_child_nodes, new_capacity, "right_child_nodes", "resize_transition_label_arrays")
            call ralloc (sorted_integral_labels, new_capacity, "sorted_integral_labels", "resize_transition_label_arrays")
            !
            label_capacity = new_capacity
        !
        case (3)
            !
            !   Deallocate the storage for the arrays
            !
            call dalloc (first_coefficient_indices, "first_coefficient_indices", "resize_transition_label_arrays")
            call dalloc (node_integral_labels, "node_integral_labels", "resize_transition_label_arrays")
            call dalloc (parent_nodes, "parent_nodes", "resize_transition_label_arrays")
            call dalloc (left_child_nodes, "left_child_nodes", "resize_transition_label_arrays")
            call dalloc (coefficient_end_offsets, "coefficient_end_offsets", "resize_transition_label_arrays")
            call dalloc (right_child_nodes, "right_child_nodes", "resize_transition_label_arrays")
            call dalloc (sorted_integral_labels, "sorted_integral_labels", "resize_transition_label_arrays")
        !
        case default
            !
            write (output_unit, *) "resize_transition_label_arrays: Invalid resize_mode = ", resize_mode
            stop
            !
        end select
        !
        return
    end subroutine resize_transition_label_arrays
end module resize_transition_label_arrays_interface
