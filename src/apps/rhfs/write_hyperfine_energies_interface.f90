module write_hyperfine_energies_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_hyperfine_energies (hamiltonian_energy_shift, e, jtot, ipar, level_indices, level_count, mode)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: hamiltonian_energy_shift
            integer, intent(in) :: level_count
            real(real64), intent(in) :: e(level_count)
            integer, intent(in) :: jtot(level_count)
            integer, intent(in) :: ipar(level_count)
            integer, intent(in) :: level_indices(level_count)
            integer, intent(in) :: mode
        end subroutine write_hyperfine_energies
    end interface
end module write_hyperfine_energies_interface
