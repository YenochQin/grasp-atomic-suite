! Select a hydrogenic self-energy factor for the requested orbital.
module hydrogenic_self_energy_factor_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: hydrogenic_self_energy_factor
contains

    !   An estimate of the function  F (Z*\alpha) is computed here.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 09 Oct 1990
    real(real64) function hydrogenic_self_energy_factor (principal_n, kappa, nuclear_charge)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use mohr_self_energy_factor_interface, only: mohr_self_energy_factor
        use klarsfeld_maquet_self_energy_factor_interface, only: klarsfeld_maquet_self_energy_factor
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: principal_n
        integer, intent(in) :: kappa
        real(real64)  :: nuclear_charge

        ! Local declarations
        integer :: effective_principal_n
        real(real64) :: unscaled_self_energy_factor

        if (principal_n <= 2) then
            call mohr_self_energy_factor (principal_n, kappa, nuclear_charge, unscaled_self_energy_factor)
        else
            if (kappa==(-1) .or. kappa==1 .or. kappa==(-2)) then
                effective_principal_n = 2
                call mohr_self_energy_factor (effective_principal_n, kappa, nuclear_charge, unscaled_self_energy_factor)
            else
                call klarsfeld_maquet_self_energy_factor (principal_n, kappa, nuclear_charge, unscaled_self_energy_factor)
            end if
        end if
        hydrogenic_self_energy_factor = unscaled_self_energy_factor
        return
    end function hydrogenic_self_energy_factor
end module hydrogenic_self_energy_factor_interface
