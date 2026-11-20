module csf_pair_subshell_couplings
    ! Module dependencies
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: bra_background_subshell_codes(3, max_orbitals), ket_background_subshell_codes(3, max_orbitals)
    integer :: bra_active_pair_codes(3), ket_active_pair_codes(3)
end module csf_pair_subshell_couplings
