module build_subshell_transformation_matrix_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_subshell_transformation_matrix (x, t, work, norb)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: norb
            real(real64) :: x(norb, norb)
            real(real64), intent(out) :: t(norb, norb)
            real(real64), intent(inout) :: work(*)
        end subroutine build_subshell_transformation_matrix
    end interface
end module build_subshell_transformation_matrix_interface
