module mohr_self_energy_factor_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine mohr_self_energy_factor (principal_n, kappa, nuclear_charge, self_energy_factor)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Local constants
            integer, parameter :: numval = 12

            ! Arguments
            integer, intent(in) :: principal_n
            integer, intent(in) :: kappa
            real(real64) :: nuclear_charge
            real(real64), intent(out) :: self_energy_factor
        end subroutine mohr_self_energy_factor
    end interface
end module mohr_self_energy_factor_interface
