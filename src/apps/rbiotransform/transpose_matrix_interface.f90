module transpose_matrix_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine transpose_matrix (xin, nrow, ncol, xout)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nrow
            integer, intent(in) :: ncol
            real(real64), intent(in) :: xin(nrow, ncol)
            real(real64), intent(out) :: xout(ncol, nrow)
        end subroutine transpose_matrix
    end interface
end module transpose_matrix_interface
