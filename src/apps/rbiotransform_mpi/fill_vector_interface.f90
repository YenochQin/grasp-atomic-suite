module fill_vector_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine fill_vector (vector, value, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ndim
            real(real64), intent(out) :: vector(ndim)
            real(real64), intent(in) :: value
        end subroutine fill_vector
    end interface
end module fill_vector_interface
