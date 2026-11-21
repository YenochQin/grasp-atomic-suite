program single_subshell_tensor
    use iso_fortran_env, only: real64
    use grasp_parameters, only: max_orbitals, csf_total_j_slot, csf_parity_slot, csf_angular_code_count
    use orbital_csf_basis, only: orbital_count, csf_count, orbital_kappa, orbital_principal_n, orbital_two_j, &
        orbital_l, csf_subshell_occupations
    use csf_coupling_data, only: csf_subshell_quantum_codes, csf_angular_codes
    use csf_pair_quantum_numbers, only: prepare_csf_pair_quantum_numbers
    use one_electron_tensor_coefficients, only: calculate_one_electron_tensor_coefficients
    use log_factorial_table, only: initialize_log_factorial_table
    implicit none (type, external)
    real(real64) :: coefficients(max_orbitals), first_coefficient
    integer :: bra_orbital, ket_orbital

    ! Two 4f(7/2)^4 states (v=2,4; J=2) and a second pair with two open shells.
    orbital_count = 2
    csf_count = 3
    orbital_kappa(:2) = [-4,-2]
    orbital_principal_n(:2) = [4,5]
    orbital_two_j(:2) = [7,3]
    orbital_l(:2) = [3,1]
    allocate(csf_subshell_occupations(max_orbitals,3), csf_subshell_quantum_codes(max_orbitals,3,3))
    allocate(csf_angular_codes(csf_angular_code_count,3))
    csf_subshell_occupations = 0
    csf_subshell_occupations(1,:) = [4,4,2]
    csf_subshell_occupations(2,:) = [0,0,2]
    csf_subshell_quantum_codes = 0
    csf_subshell_quantum_codes(1,1,:) = [2,4,2]
    csf_subshell_quantum_codes(1,3,:) = 5
    csf_subshell_quantum_codes(2,3,:) = 1
    csf_angular_codes = 0
    csf_angular_codes(:2,:) = 5
    csf_angular_codes(csf_total_j_slot,:) = 5
    csf_angular_codes(csf_parity_slot,:) = 1
    call initialize_log_factorial_table
    call calculate_one_electron_tensor_coefficients(2,1,1,2,bra_orbital,ket_orbital,coefficients)
    first_coefficient = coefficients(1)
    call prepare_csf_pair_quantum_numbers(3,3)
    call calculate_one_electron_tensor_coefficients(2,1,1,2,bra_orbital,ket_orbital,coefficients)
    if (abs(first_coefficient + 1.038051549976283_real64) > 1.0e-12_real64 .or. &
        abs(coefficients(1) - first_coefficient) > 1.0e-12_real64) error stop "single-subshell rank-2 coefficient"
end program single_subshell_tensor
