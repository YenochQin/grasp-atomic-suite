! Check the damping controls before applying orbital updates.
module check_orbital_damping_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: check_orbital_damping
contains

    !   This subroutine determines the damping factor appropriate to the
    !   present  orbital. The algorithm is taken from C Froese Fischer's
    !   program MCHF.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 08 Dec 1992
    ! Modified by C. Froese Fischer. Last update: 07 Apr 2009
    subroutine check_orbital_damping(ipr, orbital_index, previous_relative_energy_change, orbital_energy_change)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use scf_damping_factors, only: orbital_damping_factors
        use orbital_csf_basis, only: orbital_energy_parameters
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: ipr
        integer, intent(in) :: orbital_index
        real(real64), intent(inout) :: previous_relative_energy_change
        real(real64), intent(inout) :: orbital_energy_change

        ! Local declarations
        logical :: adaptive_damping_enabled
        ! The damping is adaptive (i.e., can be modified by this SUBROUTINE)
        ! if and only if ODAMP(J) .GE. 0.0
        adaptive_damping_enabled = orbital_damping_factors(orbital_index) >= 0.0_real64
        if (ipr /= orbital_index) then
            if (adaptive_damping_enabled) orbital_damping_factors(orbital_index) = &
                & 0.75_real64*orbital_damping_factors(orbital_index)
        else
            orbital_energy_change = orbital_energy_change-orbital_energy_parameters(orbital_index)
            if (adaptive_damping_enabled) then
                if (previous_relative_energy_change*orbital_energy_change > 0.0_real64) then
                    orbital_damping_factors(orbital_index) = 0.75_real64*orbital_damping_factors(orbital_index)
                else
                    orbital_damping_factors(orbital_index) = 0.25_real64+0.75_real64*orbital_damping_factors(orbital_index)
                end if
            end if
        end if
        if (ipr == orbital_index) then
            previous_relative_energy_change = orbital_energy_change
        else
            previous_relative_energy_change = orbital_energy_change-orbital_energy_parameters(orbital_index)
        end if
        ipr = orbital_index
        return
    end subroutine check_orbital_damping
end module check_orbital_damping_interface
