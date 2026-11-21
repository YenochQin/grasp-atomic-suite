module try_allocate_sort_buffers_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: try_allocate_sort_buffers
contains

    !***********************************************************************
    !                                                                      *
    subroutine try_allocate_sort_buffers(coefficient_capacity_to_probe, allocation_succeeded_flag)
        ! Module dependencies
        !                                                                      *
        !   Written by G. Gaigalas            Last revision: 27 October 2017   *
        !                                                                      *
        !                                                                      *
        !***********************************************************************
        !
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: dalloc
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: coefficient_capacity_to_probe
        integer, intent(out) :: allocation_succeeded_flag

        ! Local declarations
        real(real64), pointer :: coefficient_buffer_probe(:)
        integer, pointer :: integer_buffer_probe(:)
        integer :: allocation_status
        !-----------------------------------------------
        allocation_succeeded_flag = 1
        allocate(coefficient_buffer_probe(coefficient_capacity_to_probe),stat=allocation_status)
        if (allocation_status /= 0) then
            allocation_succeeded_flag = 0
            return
        end if
        allocate(integer_buffer_probe(3*coefficient_capacity_to_probe),stat=allocation_status)
        if (allocation_status /= 0) then
            allocation_succeeded_flag = 0
            call dalloc (coefficient_buffer_probe,"coefficient_buffer_probe","try_allocate_sort_buffers")
            return
        end if
        call dalloc (coefficient_buffer_probe,"coefficient_buffer_probe","try_allocate_sort_buffers")
        call dalloc (integer_buffer_probe,"integer_buffer_probe","try_allocate_sort_buffers")
        return
    end subroutine try_allocate_sort_buffers
end module try_allocate_sort_buffers_interface
