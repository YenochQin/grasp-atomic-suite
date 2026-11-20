!***********************************************************************
!                                                                      *
real(real64) function hydrogenic_self_energy_factor (principal_n, kappa, nuclear_charge)
    ! Module dependencies
    !                                                                      *
    !   An estimate of the function  F (Z*\alpha) is computed here.        *
    !                                                                      *
    !   Call(s) to: [RCI92]: KLAMAQ, MOHR.                                 *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1990   *
    !                                                                      *
    !***********************************************************************
    !
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
    integer :: neff
    real(real64) :: value
    !-----------------------------------------------
    !
    if (principal_n <= 2) then
        call mohr_self_energy_factor (principal_n, kappa, nuclear_charge, value)
    else
        if (kappa==(-1) .or. kappa==1 .or. kappa==(-2)) then
            neff = 2
            call mohr_self_energy_factor (neff, kappa, nuclear_charge, value)
        else
            call klarsfeld_maquet_self_energy_factor (principal_n, kappa, nuclear_charge, value)
        end if
    end if
    !
    hydrogenic_self_energy_factor = value
    !
    return
end function hydrogenic_self_energy_factor
