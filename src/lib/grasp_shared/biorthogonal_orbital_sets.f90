module biorthogonal_orbital_sets
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    character :: final_orbital_symmetry_labels(max_orbitals)*2, initial_orbital_symmetry_labels(max_orbitals)*2
    real(real64) :: final_orbital_energy_parameters(max_orbitals), final_orbital_origin_exponents(max_orbitals) , &
        initial_orbital_energy_parameters(max_orbitals), initial_orbital_origin_exponents(max_orbitals)
    integer :: final_csf_count, final_orbital_count, initial_csf_count, initial_orbital_count
    integer :: final_orbital_principal_n(max_orbitals), final_orbital_kappa(max_orbitals), &
        initial_orbital_principal_n(max_orbitals), initial_orbital_kappa(max_orbitals)
    integer :: final_orbital_l(max_orbitals), final_orbital_two_j(max_orbitals), initial_orbital_l(max_orbitals), &
        initial_orbital_two_j(max_orbitals)
    integer :: merged_to_final_orbital_index(max_orbitals), merged_to_initial_orbital_index(max_orbitals)
end module biorthogonal_orbital_sets
