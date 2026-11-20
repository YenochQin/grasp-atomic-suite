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
            call alloc (first_coefficient_indices, label_capacity, "IBEG", "ALCLLA")
            call alloc (node_integral_labels, label_capacity, "ILAB", "ALCLLA")
            call alloc (parent_nodes, label_capacity, "ILAST", "ALCLLA")
            call alloc (left_child_nodes, label_capacity, "ILEFT", "ALCLLA")
            call alloc (coefficient_end_offsets, label_capacity, "IPTCSF", "ALCLLA")
            call alloc (right_child_nodes, label_capacity, "IRIGHT", "ALCLLA")
            call alloc (sorted_integral_labels, label_capacity, "LBLINT", "ALCLLA")
        !
        case (2)
            !
            !   Double the allocation of storage for the arrays
            !
            new_capacity = 2*label_capacity
            !
            call ralloc (first_coefficient_indices, new_capacity, "IBEG", "ALCLLA")
            call ralloc (node_integral_labels, new_capacity, "ILAB", "ALCLLA")
            call ralloc (parent_nodes, new_capacity, "ILAST", "ALCLLA")
            call ralloc (left_child_nodes, new_capacity, "ILEFT", "ALCLLA")
            call ralloc (coefficient_end_offsets, new_capacity, "IPTCSF", "ALCLLA")
            call ralloc (right_child_nodes, new_capacity, "IRIGHT", "ALCLLA")
            call ralloc (sorted_integral_labels, new_capacity, "LBLINT", "ALCLLA")
            !
            label_capacity = new_capacity
        !
        case (3)
            !
            !   Deallocate the storage for the arrays
            !
            call dalloc (first_coefficient_indices, "IBEG", "ALCLLA")
            call dalloc (node_integral_labels, "ILAB", "ALCLLA")
            call dalloc (parent_nodes, "ILAST", "ALCLLA")
            call dalloc (left_child_nodes, "ILEFT", "ALCLLA")
            call dalloc (coefficient_end_offsets, "IPTCSF", "ALCLLA")
            call dalloc (right_child_nodes, "IRIGHT", "ALCLLA")
            call dalloc (sorted_integral_labels, "LBLINT", "ALCLLA")
        !
        case default
            !
            write (output_unit, *) "ALCLLA: Invalid argument IMODE = ", resize_mode
            stop
            !
        end select
        !
        return
    end subroutine resize_transition_label_arrays
end module resize_transition_label_arrays_interface
