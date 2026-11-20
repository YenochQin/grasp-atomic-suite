! Allocate, grow, or release the shared buffer for angular coefficient values and labels.
module angular_coefficient_buffer
    implicit none (type, external)

    ! Module interface
    private
    public :: alcbuf
contains
    subroutine alcbuf(buffer_operation)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit

        use array_memory_management, only: alloc, dalloc, ralloc
        use angular_coefficient_buffer_state,       only: angular_coefficient_labels, angular_coefficients, angular_buffer_capacity
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: buffer_operation

        ! Local declarations
        integer :: new_buffer_capacity, ierr
        select case (buffer_operation)
        case (1)
            angular_buffer_capacity = 10
            call alloc (angular_coefficient_labels, 6,angular_buffer_capacity, "LABEL", "ALCBUF")
            call alloc (angular_coefficients, angular_buffer_capacity, "COEFF", "ALCBUF")
        case (2)
            new_buffer_capacity = angular_buffer_capacity + angular_buffer_capacity
            call ralloc(angular_coefficient_labels,6,new_buffer_capacity, "LABEL", "ALCBUF")
            call ralloc(angular_coefficients,new_buffer_capacity, "COEFF", "ALCBUF")
            angular_buffer_capacity = new_buffer_capacity
        case (3)
            call dalloc(angular_coefficient_labels, "LABEL", "ALCBUF")
            call dalloc(angular_coefficients, "COEFF", "ALCBUF")
        case default
            write (output_unit, *) "ALCBUF: Invalid argument MODE = ", buffer_operation
            stop
        end select
        return
    end subroutine alcbuf
end module angular_coefficient_buffer
