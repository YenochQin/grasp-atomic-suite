module electron_screening_factor_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: electron_screening_factor
contains

    !***********************************************************************
    !                                                                      *
    real(real64) function electron_screening_factor (j)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real32
        !                                                                      *
        !   This  routine estimates the screening (or shielding) for orbital   *
        !   J according to the relativistic hydrogenic energy formula.         *
        !                                                                      *
        !   Written by Yu Zou, at Vanderbilt University: 02 Mar 2000           *
        !                                                                      *
        !   Formula:                                                           *
        !    epsilon = E(n,k)/Z^2_{eff}                                        *
        !    beta = alpha Z_{eff}                                              *
        !    epsilon = 1/beta^{2} [(1+beta^2/nu^2)^{-1/2} - 1]                 *
        !    nu = n + (k^2 - beta^2)^{1/2} - |k|                               *
        !   where E(n,k) is the orbital energy for the specific principal      *
        !   quantum n and k quantum number, alpha is the fine structure        *
        !   constant.                                                          *
        !   Note that SHIELD is set to -C if Zeff is out of range. This will   *
        !   ensure that the self-energy for such orbital is zero.              *
        !   Reference: M J Seaton, Rep. Prog. Phys., Vol.46, P167,1983         *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: nuclear_charge, speed_of_light_au
        use orbital_csf_basis, only: orbital_energy_parameters, orbital_kappa, orbital_principal_n
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: j

        ! Local declarations
        real(real64) :: a, aa, bb, cc, t, beta
        !-----------------------------------------------
        !
        !
        ! a=[(1+\alpha^2 E)^{-2} - 1]^{-1/2}
        ! note E(J) = -E(n,k)
        a = 1 - orbital_energy_parameters(j)/(speed_of_light_au*speed_of_light_au)
        a = 1/(a*a) - 1
        if (a < 0.0_real32) then
            electron_screening_factor = -speed_of_light_au
            return
        end if
        a = 1/sqrt(a)
        !
        aa = 1 + a*a
        bb = a*dble(orbital_principal_n(j)-iabs(orbital_kappa(j)))
        cc = dble((orbital_principal_n(j)-iabs(orbital_kappa(j)))**2-orbital_kappa(j)**2)
        t = bb*bb - aa*cc
        if (t < 0.0_real32) then
            electron_screening_factor = -speed_of_light_au
            return
        end if
        beta = bb + sqrt(t)
        beta = beta/aa
        electron_screening_factor = nuclear_charge - beta*speed_of_light_au
        !     print *, j, '  Initial Screen Factor = ',shield
        return
    end function electron_screening_factor
end module electron_screening_factor_interface
