! Holds occupations, coupling data, and the peel/core orbital lists for the current CSF pair.
module csf_pair_angular_state
    ! Module dependencies
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: bra_occupations(max_orbitals), ket_occupations(max_orbitals)
    ! Each intermediate coupling chain is in peel order, with a node for each shell after the first.
    integer :: bra_intermediate_couplings(max_orbitals), ket_intermediate_couplings(max_orbitals)
    ! These pair-local arrays put the quantum-number selector first for each orbital.
    integer :: bra_subshell_quantum_codes(3, max_orbitals), ket_subshell_quantum_codes(3, max_orbitals)
    integer :: peel_orbital_indices(max_orbitals), core_orbital_indices(max_orbitals)
    integer :: peel_orbital_count, core_orbital_count
end module csf_pair_angular_state
