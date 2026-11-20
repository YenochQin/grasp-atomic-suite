module combine_vectors_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine combine_vectors (c, a, b, faca, facb, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ndim
            real(real64), intent(in)  :: a(ndim)
            real(real64), intent(in)  :: b(ndim)
            real(real64), intent(out) :: c(ndim)
            real(real64), intent(in) :: faca
            real(real64), intent(in) :: facb
        end subroutine combine_vectors
    end interface
end module combine_vectors_interface
