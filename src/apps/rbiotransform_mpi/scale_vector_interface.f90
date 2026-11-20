module scale_vector_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine scale_vector (vector, factor, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ndim
            real(real64), intent(inout) :: vector(ndim)
            real(real64), intent(in) :: factor
        end subroutine scale_vector
    end interface
end module scale_vector_interface
