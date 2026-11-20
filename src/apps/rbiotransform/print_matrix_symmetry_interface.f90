module print_matrix_symmetry_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine print_matrix_symmetry (a, matdim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: a(1)
            integer, intent(in) :: matdim
            !...This routine performs I/O.
        end subroutine print_matrix_symmetry
    end interface
end module print_matrix_symmetry_interface
