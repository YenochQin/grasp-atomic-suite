module scale_vector_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: scale_vector
contains

    ! CALCULATE SCALAR(FACTOR) TIMES VECTOR
    subroutine scale_vector(vector, factor, ndim)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: ndim
        real(real64), intent(in) :: factor
        real(real64), intent(inout) :: vector(ndim)

        vector(:ndim) = vector(:ndim)*factor
        return
    end subroutine scale_vector
end module scale_vector_interface
