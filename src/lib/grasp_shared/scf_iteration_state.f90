module scf_iteration_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,    only: max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: orbital_occupation_weights(max_orbitals)
    integer :: radial_solution_methods(max_orbitals)
    real(real64) :: orbital_self_consistency_errors(max_orbitals)
    real(real64) :: lower_trial_energy_bound, upper_trial_energy_bound, initial_lower_trial_energy_bound, &
        initial_upper_trial_energy_bound, asymptotic_nuclear_charge
    integer :: direct_coefficient_count, exchange_coefficient_count, y_potential_coefficient_count, &
        direct_coefficient_capacity, exchange_coefficient_capacity, y_potential_coefficient_capacity
    real(real64), pointer :: direct_coefficients(:), exchange_coefficients(:), y_potential_coefficients(:)
    integer, pointer :: direct_orbital_indices(:), exchange_coefficient_labels(:), y_potential_coefficient_labels(:)
end module scf_iteration_state
