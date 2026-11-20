module write_scf_energies_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_scf_energies(e,level_indices,level_count,mode)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: level_count
            real(real64), intent(in) :: e(level_count)
            integer, intent(in) :: level_indices(level_count)
            integer, intent(in) :: mode
        end subroutine write_scf_energies
    end interface
end module write_scf_energies_interface
