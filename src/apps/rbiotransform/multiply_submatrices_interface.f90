module multiply_submatrices_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine multiply_submatrices (c, a, b, ncrow, nccol, narow, nacol, nbrow, nbcol&
            , itrnsp)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ncrow
            integer, intent(in) :: nccol
            real(real64), intent(inout) :: c(ncrow, nccol)
            integer, intent(in) :: narow
            integer, intent(in) :: nacol
            real(real64), intent(in) :: a(narow, nacol)
            integer, intent(in) :: nbrow
            integer, intent(in) :: nbcol
            real(real64), intent(in) :: b(nbrow, nbcol)
            integer, intent(in) :: itrnsp
        end subroutine multiply_submatrices
    end interface
end module multiply_submatrices_interface
