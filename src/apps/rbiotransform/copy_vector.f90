module copy_vector_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: copy_vector
contains

    subroutine copy_vector(from, to, ndim)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: ndim
        real(real64), intent(in)  :: from(ndim)
        real(real64), intent(out) :: to(ndim)
        to(:ndim) = from(:ndim)
        return
    end subroutine copy_vector
end module copy_vector_interface
