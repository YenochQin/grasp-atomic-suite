module check_orbital_damping_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: check_orbital_damping
contains

    !***********************************************************************
    !                                                                      *
    subroutine check_orbital_damping(ipr, orbital_index, previous_relative_energy_change, orbital_energy_change)
        ! Module dependencies
        !                                                                      *
        !   This subroutine determines the damping factor appropriate to the   *
        !   present  orbital. The algorithm is taken from C Froese Fischer's   *
        !   program MCHF.                                                      *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use scf_damping_factors, only: orbital_damping_factors
        use orbital_csf_basis, only: orbital_energy_parameters, previous_orbital_energy_changes
        implicit none (type, external)

        ! Arguments
        integer, intent(inout) :: ipr
        integer, intent(in) :: orbital_index
        real(real64), intent(inout) :: previous_relative_energy_change
        real(real64), intent(inout) :: orbital_energy_change

        ! Local declarations
        logical :: adaptive_damping_enabled
        !-----------------------------------------------
        !
        !   The damping is adaptive (i.e., can be modified by this SUBROUTINE)
        !   if and only if ODAMP(J) .GE. 0.0
        !
        adaptive_damping_enabled = orbital_damping_factors(orbital_index) >= 0.0_real64
        !
        orbital_energy_change = (orbital_energy_change-orbital_energy_parameters(orbital_index))/orbital_energy_change
        if (adaptive_damping_enabled) then
            if (previous_relative_energy_change*orbital_energy_change < -1.0e-4_real64) then
                !           We have a significant oscillation, damp
                orbital_damping_factors(orbital_index) = 0.10_real64 + 0.90_real64*orbital_damping_factors(orbital_index)
            !           Print *, 'Increasing ODAMP', odamp(j)
            else
                orbital_damping_factors(orbital_index) = 0.50_real64*orbital_damping_factors(orbital_index)
            end if
        end if
        !
        !     Save the relative difference for the next update of this orbital
        previous_orbital_energy_changes(orbital_index) = orbital_energy_change
        !
        return
    end subroutine check_orbital_damping
end module check_orbital_damping_interface
