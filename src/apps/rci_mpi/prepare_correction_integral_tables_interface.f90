module prepare_correction_integral_tables_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine prepare_correction_integral_tables (max_orbital_two_j, atwinv)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: max_orbital_two_j
            real(real64), intent(out) :: atwinv
        end subroutine prepare_correction_integral_tables
    end interface
end module prepare_correction_integral_tables_interface
