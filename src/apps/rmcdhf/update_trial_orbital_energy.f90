module update_trial_orbital_energy_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: update_trial_orbital_energy
contains

    !   This  subroutine implements Part 2 of Algorithm 7.1 in  C Froese
    !   Fischer,  Comput  Phys  Rep, 3 (1986) 273-326. (The present code
    !   actually  implements  the version used in her program MCHF where
    !   differences occur.)
    !
    ! Written by Farid A Parpia, at Oxford. Last Update: 08 Dec 1992
    subroutine update_trial_orbital_energy(orbital_index, first_oscillation_sign, effective_principal_n, &
        & node_count_difference, energy_correction, &
        & solve_failed, orbital_sign_reversal_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_dirac_workspace, only: trial_orbital_last_grid_index, trial_large_component, trial_large_origin_parameter, &
            trial_small_component
        use orbital_csf_basis, only: orbital_energy_parameters, orbital_principal_n
        use scf_iteration_state, only: lower_trial_energy_bound, upper_trial_energy_bound, initial_lower_trial_energy_bound, &
            initial_upper_trial_energy_bound

        ! Procedure interfaces
        use energy_outside_bounds_interface, only: energy_outside_bounds
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index
        integer , intent(in) :: effective_principal_n
        integer , intent(in) :: node_count_difference
        integer , intent(out) :: orbital_sign_reversal_count
        real(real64) , intent(in) :: first_oscillation_sign
        real(real64) , intent(inout) :: energy_correction
        logical , intent(out) :: solve_failed

        ! Local constants
        real(real64), parameter :: single_node_energy_change_threshold = 2.0e-2_real64
        real(real64), parameter :: near_energy_bound_threshold = 5.0e-2_real64
        real(real64), parameter :: minimum_absolute_energy_correction = 1.0e-3_real64
        real(real64), parameter :: minimum_relative_energy_correction = 1.0e-5_real64
        real(real64), parameter :: node_energy_scaling_exponent = 2.5_real64

        ! Local declarations
        real(real64) :: current_orbital_energy, absolute_energy_correction, relative_energy_correction, principal_n, &
            & trial_principal_n, node_energy_scale, trial_energy, relative_bound_distance
        ! Determine if the iterative process has succeeded
        if (first_oscillation_sign>0.0_real64 .and. node_count_difference==0) then
            solve_failed = .false.
        else
            solve_failed = .true.
        end if
        ! Inversion counter
        orbital_sign_reversal_count = 0
        ! If unsuccessful, obtain a new estimate of the eigenvalue
        if (solve_failed) then
            ! Define quantities used later
            current_orbital_energy = orbital_energy_parameters(orbital_index)
            absolute_energy_correction = abs(energy_correction)
            relative_energy_correction = abs(energy_correction/current_orbital_energy)
            principal_n = dble(orbital_principal_n(orbital_index))
            trial_principal_n = dble(effective_principal_n)
            node_energy_scale = (trial_principal_n/principal_n)**node_energy_scaling_exponent
            if (abs(node_count_difference)==1 .and. relative_energy_correction>single_node_energy_change_threshold .or. &
                & node_count_difference==0 .and. relative_energy_correction>=minimum_relative_energy_correction .and. &
                absolute_energy_correction>=minimum_absolute_energy_correction) then
                do
                    trial_energy = current_orbital_energy + energy_correction
                    if (energy_outside_bounds(trial_energy) .and. node_count_difference/=0) trial_energy = &
                        & current_orbital_energy*node_energy_scale
                    if (energy_outside_bounds(trial_energy)) trial_energy = current_orbital_energy - energy_correction
                    if (energy_outside_bounds(trial_energy)) then
                        energy_correction = 0.5_real64*energy_correction
                    else
                        exit
                    end if
                end do
            else if (node_count_difference == 0) then
                trial_energy = current_orbital_energy
                orbital_sign_reversal_count = 1
                trial_large_origin_parameter = -trial_large_origin_parameter
                trial_large_component(:trial_orbital_last_grid_index) = -trial_large_component(:trial_orbital_last_grid_index)
                trial_small_component(:trial_orbital_last_grid_index) = -trial_small_component(:trial_orbital_last_grid_index)
                solve_failed = .false.
            else if (first_oscillation_sign < 0.0_real64) then
                trial_energy = current_orbital_energy*node_energy_scale
                if (energy_outside_bounds(trial_energy)) trial_energy = 0.5_real64*(upper_trial_energy_bound + &
                    current_orbital_energy)
                if (energy_outside_bounds(trial_energy)) trial_energy = 0.5_real64*(lower_trial_energy_bound + &
                    current_orbital_energy)
            else if (node_count_difference < 0) then
                relative_bound_distance = 1.0_real64 - current_orbital_energy/upper_trial_energy_bound
                upper_trial_energy_bound = current_orbital_energy
                if (relative_bound_distance < near_energy_bound_threshold) then
                    initial_upper_trial_energy_bound = initial_upper_trial_energy_bound*node_energy_scale
                else
                    initial_upper_trial_energy_bound = current_orbital_energy*node_energy_scale
                end if
                trial_energy = initial_upper_trial_energy_bound
                if (energy_outside_bounds(trial_energy)) trial_energy = 0.5_real64*(upper_trial_energy_bound + &
                    lower_trial_energy_bound)
            else if (node_count_difference > 0) then
                relative_bound_distance = 1.0_real64 - lower_trial_energy_bound/current_orbital_energy
                lower_trial_energy_bound = current_orbital_energy
                if (relative_bound_distance < near_energy_bound_threshold) then
                    initial_lower_trial_energy_bound = initial_lower_trial_energy_bound*node_energy_scale
                else
                    initial_lower_trial_energy_bound = current_orbital_energy*node_energy_scale
                end if
                trial_energy = initial_lower_trial_energy_bound
                if (energy_outside_bounds(trial_energy)) trial_energy = 0.5_real64*(upper_trial_energy_bound + &
                    lower_trial_energy_bound)
            end if
            orbital_energy_parameters(orbital_index) = trial_energy
        end if
        return
    end subroutine update_trial_orbital_energy
end module update_trial_orbital_energy_interface
