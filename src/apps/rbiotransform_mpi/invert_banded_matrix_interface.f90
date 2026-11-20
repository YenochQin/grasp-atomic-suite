module invert_banded_matrix_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine invert_banded_matrix (a, el, n, determ, epsil, itest, nsize)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nsize
            real(real64), intent(inout) :: a(nsize, nsize)
            real(real64), intent(inout) :: el(nsize, nsize)
            integer, intent(in) :: n
            real(real64), intent(out) :: determ
            real(real64), intent(in) :: epsil
            integer, intent(out) :: itest
        end subroutine invert_banded_matrix
    end interface
end module invert_banded_matrix_interface
