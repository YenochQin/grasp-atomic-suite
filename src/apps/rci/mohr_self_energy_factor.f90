! Evaluate the Mohr hydrogenic self-energy factor.
module mohr_self_energy_factor_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: mohr_self_energy_factor
contains

    !   The  function  F (Z*alpha)  for the  1s  2s  2p-  2p  symmetries
    !   is computed here.    A value is obtained by interpolating in, or
    !   extrapolating from, the table due to  P J Mohr.   See  P J Mohr,
    !   At Data Nucl Data Tables 29 (1983) 453.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 09 Oct 1992
    subroutine mohr_self_energy_factor(principal_n, kappa, nuclear_charge, self_energy_factor)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use aitken_interpolation, only: interpolate_aitken
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: principal_n
        integer , intent(in) :: kappa
        real(real64)  :: nuclear_charge
        real(real64) , intent(out) :: self_energy_factor

        ! Local constants
        !   Number of data points
        integer, parameter :: tabulated_charge_count = 12

        ! Local declarations
        real(real64) :: factors_1s(tabulated_charge_count), factors_2s(tabulated_charge_count), &
            & factors_2p_half(tabulated_charge_count), factors_2p_three_halves(tabulated_charge_count), &
            & tabulated_nuclear_charges(tabulated_charge_count)
        real(real64) :: interpolation_tolerance, interpolated_factor

        ! Initial values
        !
        !   1s data :
        data factors_1s/ 10.3168_real64, 4.6540_real64, 3.2460_real64, 2.5519_real64, 2.1351_real64, &
            1.8644_real64, 1.6838_real64, 1.5675_real64, 1.5032_real64, 1.4880_real64, 1.5317_real64, &
            1.6614_real64/
        ! 2s data:
        data factors_2s/ 10.5468_real64, 4.8930_real64, 3.5063_real64, 2.8391_real64, 2.4550_real64, &
            2.2244_real64, 2.0948_real64, 2.0435_real64, 2.0650_real64, 2.1690_real64, 2.3870_real64, &
            2.7980_real64/
        ! 2p- data:
        data factors_2p_half/ -0.1264_real64, -0.1145_real64, -0.0922_real64, -0.0641_real64, -0.0308_real64, &
            0.0082_real64, 0.0549_real64, 0.1129_real64, 0.1884_real64, 0.2934_real64, 0.4530_real64, &
            0.7250_real64/
        ! 2p data:
        data factors_2p_three_halves/ 0.1235_real64, 0.1303_real64, 0.1436_real64, 0.1604_real64, 0.1794_real64, &
            0.1999_real64, 0.2215_real64, 0.2440_real64, 0.2671_real64, 0.2906_real64, 0.3141_real64, &
            0.3367_real64/
        ! Z data:
        data tabulated_nuclear_charges/ 1.0_real64, 10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64, 60.0_real64, &
            70.0_real64, 80.0_real64, 90.0_real64, 100.0_real64, 110.0_real64/
        ! Convergence criterion for interpolation
        data interpolation_tolerance/ 1.0e-3_real64/
        ! Interpolate or issue error message as appropriate
        if (principal_n == 1) then
            if (kappa == (-1)) then
                call interpolate_aitken (tabulated_nuclear_charges, factors_1s, tabulated_charge_count, nuclear_charge, &
                    & interpolated_factor, interpolation_tolerance)
            else
                write (*, '("mohr_self_energy_factor:")')
                write (*, '(" Principal quantum number, ",i12,", kappa, ",1i3,".")') principal_n, kappa
                stop
            end if
        else if (principal_n == 2) then
            select case (kappa)
            case (-1)
                call interpolate_aitken (tabulated_nuclear_charges, factors_2s, tabulated_charge_count, nuclear_charge, &
                    & interpolated_factor, interpolation_tolerance)
            case (1)
                call interpolate_aitken (tabulated_nuclear_charges, factors_2p_half, tabulated_charge_count, nuclear_charge, &
                    & interpolated_factor, interpolation_tolerance)
            case (-2)
                call interpolate_aitken (tabulated_nuclear_charges, factors_2p_three_halves, tabulated_charge_count, &
                    & nuclear_charge, interpolated_factor, interpolation_tolerance)
            case default
                write (*, '("mohr_self_energy_factor:")')
                write (*, '(" Principal quantum number, ",i12,", kappa, ",1i3,".")') principal_n, kappa
                stop
            end select
        else
            write (*, '("mohr_self_energy_factor:")')
            write (*, '(" Principal quantum number, ",1i2,", Should be either 1 or 2.")') principal_n
            stop
        end if
        self_energy_factor = interpolated_factor
        return
        return
    end subroutine mohr_self_energy_factor
end module mohr_self_energy_factor_interface
