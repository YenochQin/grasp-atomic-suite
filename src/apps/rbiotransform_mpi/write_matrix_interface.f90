module write_matrix_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_matrix (a, nrow, ncol, nmrow, nmcol)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nmrow
            integer, intent(in) :: nmcol
            real(real64), intent(in) :: a(nmrow, nmcol)
            integer, intent(in) :: nrow
            integer, intent(in) :: ncol
            !...This routine performs I/O.
        end subroutine write_matrix
    end interface
end module write_matrix_interface
