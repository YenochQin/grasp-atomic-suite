module copy_vector_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine copy_vector (from, to, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ndim
            real(real64), intent(in) :: from(ndim)
            real(real64), intent(out) :: to(ndim)
        end subroutine copy_vector
    end interface
end module copy_vector_interface
