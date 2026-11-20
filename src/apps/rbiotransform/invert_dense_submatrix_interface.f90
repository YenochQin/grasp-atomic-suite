module invert_dense_submatrix_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine invert_dense_submatrix (a, b, matdim, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: matdim
            real(real64), intent(inout) :: a(matdim, matdim)
            real(real64) :: b(matdim, matdim)
            integer, intent(in) :: ndim
        end subroutine invert_dense_submatrix
    end interface
end module invert_dense_submatrix_interface
