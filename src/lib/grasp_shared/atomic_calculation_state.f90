! Shares atomic-system parameters, unit conversions, SCF controls, and correction arrays.
module atomic_calculation_state
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only:  max_radial_points
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64) :: maximum_decimal_exponent, maximum_exp_argument, minimum_exp_argument, machine_precision
    real(real64) :: hartree_to_wavenumber_cm, hartree_to_ev, speed_of_light_cm_per_s, einstein_a_au_to_si, einstein_b_au_to_si
    real(real64) :: femtometer_to_bohr, electron_mass_u
    integer :: ionic_charge, electron_count
    real(real64) :: nuclear_mass_in_electron_masses, nuclear_charge
    integer :: final_ionic_charge, final_electron_count
    real(real64) :: final_nuclear_mass_in_electron_masses, final_nuclear_charge
    integer :: initial_ionic_charge, initial_electron_count
    real(real64) :: initial_nuclear_mass_in_electron_masses, initial_nuclear_charge
    integer :: left_electron_count
    real(real64) :: speed_of_light_au
    real(real64) :: proton_mass_u, bohr_radius_cm
    integer :: maximum_scf_iterations, orbital_updates_per_scf_iteration, minimum_radial_solver_iterations
    real(real64) :: working_tolerance
    real(real64), pointer :: level_average_weights(:), requested_level_weights(:)
    integer :: selected_level_count, highest_requested_level_index
    integer, pointer :: requested_level_indices(:)
    real(real64) :: vacuum_speed_of_light_au, pi
    real(real64) :: large_deferred_correction(max_radial_points), small_deferred_correction(max_radial_points)
end module atomic_calculation_state
