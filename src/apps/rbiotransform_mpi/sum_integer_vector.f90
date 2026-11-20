module sum_integer_vector_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: sum_integer_vector
contains

    !
    !     ------------------------------------------------------------------
    !     I E L S U M
    !     ------------------------------------------------------------------
    !
    integer function sum_integer_vector (ivec, nelmnt)
        !
        ! Sum elements of integer array
        !
        !************************************************************************
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nelmnt
        integer, intent(in) :: ivec(*)

        ! Local declarations
        integer :: isum
        !-----------------------------------------------
        !
        isum = 0
        isum = sum(ivec(:nelmnt))
        !
        sum_integer_vector = isum
        !
        return
    end function sum_integer_vector
end module sum_integer_vector_interface
