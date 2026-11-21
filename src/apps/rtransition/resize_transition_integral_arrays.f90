module resize_transition_integral_arrays_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: resize_transition_integral_arrays
contains

    !***********************************************************************
    !                                                                      *
    subroutine resize_transition_integral_arrays(transition_bra_csf_indices, transition_ket_csf_indices, &
        & transition_angular_coefficients, &
        & integral_capacity, resize_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   This  subprogram allocates (IMODE = 1), reallocates (IMODE = 2),   *
        !   and deallocates (IMODE = 3) storage for  certain arrays that are   *
        !   used in OSCL.                                                      *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, DALLOC, RALLOC.                        *
        !                                                                      *
        !   Farid A. Parpia.                      Last revision: 28 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use array_memory_management, only: alloc, dalloc, ralloc
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: integral_capacity
        integer, intent(in) :: resize_mode
        integer, pointer :: transition_bra_csf_indices(:), transition_ket_csf_indices(:)
        real(real64), pointer :: transition_angular_coefficients(:)

        ! Local declarations
        integer :: new_capacity
        !-----------------------------------------------
        !
        select case (resize_mode)
        case (1)
            !
            !   Initial array dimension
            !
            integral_capacity = 1
            !
            !   Allocate storage for arrays
            !
            call alloc (transition_bra_csf_indices, integral_capacity, "transition_bra_csf_indices", &
                "resize_transition_integral_arrays" )
            call alloc (transition_ket_csf_indices, integral_capacity,"transition_ket_csf_indices", &
                "resize_transition_integral_arrays" )
            call alloc (transition_angular_coefficients, integral_capacity, "transition_angular_coefficients", &
                "resize_transition_integral_arrays" )
        !
        case (2)
            !
            !   Double the allocation of storage for the arrays
            !
            new_capacity = 2*integral_capacity
            !
            call ralloc (transition_bra_csf_indices, new_capacity, "transition_bra_csf_indices", &
                "resize_transition_integral_arrays" )
            call ralloc (transition_ket_csf_indices, new_capacity, "transition_ket_csf_indices", &
                "resize_transition_integral_arrays" )
            call ralloc (transition_angular_coefficients, new_capacity, "transition_angular_coefficients", &
                "resize_transition_integral_arrays" )
            !
            integral_capacity = new_capacity
        !
        case (3)
            !
            !   Deallocate the storage for the arrays
            !
            call dalloc (transition_bra_csf_indices, "transition_bra_csf_indices", "resize_transition_integral_arrays" )
            call dalloc (transition_ket_csf_indices, "transition_ket_csf_indices", "resize_transition_integral_arrays" )
            call dalloc (transition_angular_coefficients,  "transition_angular_coefficients", "resize_transition_integral_arrays" )
        !
        case default
            !
            write (output_unit, *) "resize_transition_integral_arrays: Invalid resize_mode = ", resize_mode
            stop
            !
        end select
        !
        return
    end subroutine resize_transition_integral_arrays
end module resize_transition_integral_arrays_interface
