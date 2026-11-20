module estimate_orbital_energy_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_orbital_energy
contains

    !***********************************************************************
    !                                                                      *
    subroutine estimate_orbital_energy(orbital_index)
        ! Module dependencies
        !                                                                      *
        !   This  subprogram implements  Part 1 of Algorithm 7.1 of C Froese   *
        !   Fischer, Comput Phys Rep, 3 (1986) 320-321.                        *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
        !                                                                      *
        !***********************************************************************
        use atomic_calculation_state, only: electron_count, nuclear_charge, speed_of_light_au
        use orbital_csf_basis, only: orbital_kappa, orbital_principal_n
        use scf_iteration_state, only: lower_trial_energy_bound, upper_trial_energy_bound, initial_lower_trial_energy_bound, &
            initial_upper_trial_energy_bound, asymptotic_nuclear_charge
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index

        ! Local declarations
        integer :: principal_n, absolute_kappa
        real(real64) :: fine_structure_constant, speed_of_light_squared, radial_quantum_number, absolute_kappa_real, &
            & kappa_squared, charge_times_fine_structure_constant, gamma, dirac_energy_over_rest_energy
        !-----------------------------------------------
        !
        !   Initializations
        !
        fine_structure_constant = 1.0_real64/speed_of_light_au
        speed_of_light_squared = speed_of_light_au*speed_of_light_au
        principal_n = orbital_principal_n(orbital_index)
        absolute_kappa = abs(orbital_kappa(orbital_index))
        radial_quantum_number = dble(principal_n - absolute_kappa)
        absolute_kappa_real = dble(absolute_kappa)
        kappa_squared = absolute_kappa_real*absolute_kappa_real
        !
        !   Set ZINF, the asymptotic charge seen by the electron
        !
        !     ZINF = DBLE (IONCTY+1)
        !
        !   Changed on 07/06/93 by WPW
        !
        asymptotic_nuclear_charge = nuclear_charge + dble((-electron_count) + 1)
        !
        !   Set the lower bound
        !
        charge_times_fine_structure_constant = asymptotic_nuclear_charge*fine_structure_constant
        if (charge_times_fine_structure_constant < absolute_kappa_real) then
            gamma = sqrt(kappa_squared - charge_times_fine_structure_constant*charge_times_fine_structure_constant)
            dirac_energy_over_rest_energy = 1.0_real64/sqrt(1.0_real64 + (charge_times_fine_structure_constant/(gamma + &
                & radial_quantum_number + 0.5_real64))**2)
            lower_trial_energy_bound = (1.0_real64 - dirac_energy_over_rest_energy)*speed_of_light_squared
        else
            lower_trial_energy_bound = 0.25_real64*speed_of_light_squared/dble(principal_n*principal_n)
        end if
        initial_lower_trial_energy_bound = lower_trial_energy_bound
        !
        !   Set the upper bound
        !
        charge_times_fine_structure_constant = nuclear_charge*fine_structure_constant
        if (charge_times_fine_structure_constant < absolute_kappa_real) then
            gamma = sqrt(kappa_squared - charge_times_fine_structure_constant*charge_times_fine_structure_constant)
            dirac_energy_over_rest_energy = 1.0_real64/sqrt(1.0_real64 + (charge_times_fine_structure_constant/(gamma + &
                & radial_quantum_number - 0.5_real64))**2)
            upper_trial_energy_bound = (1.0_real64 - dirac_energy_over_rest_energy)*speed_of_light_squared
        else
            upper_trial_energy_bound = speed_of_light_squared + speed_of_light_squared
        end if
        initial_upper_trial_energy_bound = upper_trial_energy_bound
        !
        return
    end subroutine estimate_orbital_energy
end module estimate_orbital_energy_interface
