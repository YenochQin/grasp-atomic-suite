module sum_integer_vector_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: sum_integer_vector
contains

    ! Sum elements of integer array
    integer function sum_integer_vector (ivec, nelmnt)
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nelmnt
        integer, intent(in) :: ivec(*)

        ! Local declarations
        integer :: isum

        isum = 0
        isum = sum(ivec(:nelmnt))
        sum_integer_vector = isum
        return
    end function sum_integer_vector
end module sum_integer_vector_interface
