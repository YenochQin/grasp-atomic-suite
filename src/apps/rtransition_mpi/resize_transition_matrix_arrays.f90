module resize_transition_matrix_arrays_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: resize_transition_matrix_arrays
contains

    ! This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),
    ! and deallocates (IMODE = 3) storage for  certain arrays that are
    ! used in 12_oscl.
    !
    ! Farid A. Parpia.. Last revision: 28 Dec 1992
    subroutine resize_transition_matrix_arrays(jja, jjb, hb1, hb2, hc1, hc2, hm1, &
        hm2, lab, nptr, matrix_element_capacity, resize_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use array_memory_management, only: alloc, dalloc, ralloc
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer  :: matrix_element_capacity
        integer, intent(in) :: resize_mode
        integer, pointer :: jja(:), jjb(:), lab(:), nptr(:)
        real(real64), pointer :: hb1(:), hb2(:), hc1(:), hc2(:), hm1(:), hm2(:)

        ! Local declarations
        integer :: new_capacity

        select case (resize_mode)
        case (1)
            ! Initial array dimension
            matrix_element_capacity = 1
            ! Allocate storage for arrays
            call alloc (jja, matrix_element_capacity, "JJA", "resize_transition_matrix_arrays")
            call alloc (jjb, matrix_element_capacity, "JJB", "resize_transition_matrix_arrays" )
            call alloc (hb1, matrix_element_capacity, "HB1", "resize_transition_matrix_arrays")
            call alloc (hb2, matrix_element_capacity, "HB2", "resize_transition_matrix_arrays")
            call alloc (hc1, matrix_element_capacity, "HC1", "resize_transition_matrix_arrays")
            call alloc (hc2, matrix_element_capacity, "HC2", "resize_transition_matrix_arrays")
            call alloc (hm1, matrix_element_capacity, "HM1", "resize_transition_matrix_arrays")
            call alloc (hm2, matrix_element_capacity, "HM2", "resize_transition_matrix_arrays")
            call alloc (lab, matrix_element_capacity, "LAB", "resize_transition_matrix_arrays" )
            call alloc (nptr, matrix_element_capacity, "NPTR", "resize_transition_matrix_arrays" )
        case (2)
            ! Double the allocation of storage for the arrays
            new_capacity = 2*matrix_element_capacity
            call ralloc (jja, new_capacity, "JJA", "resize_transition_matrix_arrays")
            call ralloc (jjb, new_capacity, "JJB", "resize_transition_matrix_arrays" )
            call ralloc (hb1, new_capacity, "HB1", "resize_transition_matrix_arrays")
            call ralloc (hb2, new_capacity, "HB2", "resize_transition_matrix_arrays")
            call ralloc (hc1, new_capacity, "HC1", "resize_transition_matrix_arrays")
            call ralloc (hc2, new_capacity, "HC2", "resize_transition_matrix_arrays")
            call ralloc (hm1, new_capacity, "HM1", "resize_transition_matrix_arrays")
            call ralloc (hm2, new_capacity, "HM2", "resize_transition_matrix_arrays")
            call ralloc (lab, new_capacity, "LAB", "resize_transition_matrix_arrays" )
            call ralloc (nptr, new_capacity, "NPTR", "resize_transition_matrix_arrays" )
            matrix_element_capacity = new_capacity
        case (3)
            ! Deallocate the storage for the arrays
            call dalloc (jja, "JJA", "resize_transition_matrix_arrays")
            call dalloc (jjb, "JJB", "resize_transition_matrix_arrays" )
            call dalloc (hb1, "HB1", "resize_transition_matrix_arrays")
            call dalloc (hb2, "HB2", "resize_transition_matrix_arrays")
            call dalloc (hc1, "HC1", "resize_transition_matrix_arrays")
            call dalloc (hc2, "HC2", "resize_transition_matrix_arrays")
            call dalloc (hm1, "HM1", "resize_transition_matrix_arrays")
            call dalloc (hm2, "HM2", "resize_transition_matrix_arrays")
            call dalloc (lab, "LAB", "resize_transition_matrix_arrays" )
            call dalloc (nptr, "NPTR", "resize_transition_matrix_arrays" )
        case default
            write (output_unit, *) "resize_transition_matrix_arrays: Invalid resize_mode = ", resize_mode
            stop
        end select
        return
    end subroutine resize_transition_matrix_arrays
end module resize_transition_matrix_arrays_interface
