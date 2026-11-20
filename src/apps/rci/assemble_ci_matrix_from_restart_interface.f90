module assemble_ci_matrix_from_restart_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine assemble_ci_matrix_from_restart (irestart, nelmnt_a, elsto)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64, int64
            implicit none (type, external)

            ! Arguments
            integer, intent(in)        :: irestart
            integer(int64), intent(out) :: nelmnt_a
            real(real64), intent(in)   :: elsto
        end subroutine assemble_ci_matrix_from_restart
    end interface
end module assemble_ci_matrix_from_restart_interface
