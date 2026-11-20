module write_lande_g_factors_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_lande_g_factors (name, nci, gjc_diag, dgjc_diag)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use selected_level_indices, only: nvec
            implicit none (type, external)

            ! Arguments
            character(len=24), intent(in) :: name
            integer, intent(in) :: nci
            real(real64), intent(in) :: gjc_diag(nvec), dgjc_diag(nvec)
        end subroutine write_lande_g_factors
    end interface
end module write_lande_g_factors_interface
