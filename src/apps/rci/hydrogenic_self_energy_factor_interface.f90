module hydrogenic_self_energy_factor_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function hydrogenic_self_energy_factor (principal_n, kappa, nuclear_charge)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: principal_n
            integer, intent(in) :: kappa
            real(real64) :: nuclear_charge
        end function hydrogenic_self_energy_factor
    end interface
end module hydrogenic_self_energy_factor_interface
