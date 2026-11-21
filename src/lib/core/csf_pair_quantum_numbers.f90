! Extract the occupations and coupling chains of a CSF pair into the angular workspace.
module csf_pair_quantum_numbers
    implicit none (type, external)

    ! Module interface
    private
    public :: prepare_csf_pair_quantum_numbers
contains
    subroutine prepare_csf_pair_quantum_numbers(bra_csf_index, ket_csf_index)
        ! Module dependencies
        use csf_pair_angular_state, only: bra_intermediate_couplings, ket_intermediate_couplings, bra_subshell_quantum_codes, &
            ket_subshell_quantum_codes, peel_orbital_indices, core_orbital_indices, core_orbital_count, peel_orbital_count, &
            bra_occupations, ket_occupations
        use orbital_csf_basis, only: orbital_count
        use csf_quantum_queries, only: csf_occupation, csf_subshell_quantum_code, csf_cumulative_coupling
        use subshell_occupation_class, only: classify_subshell_occupation
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index
        integer, intent(in) :: ket_csf_index

        ! Local declarations
        integer :: orbital_or_peel_index, quantum_number_selector, peel_orbital
        do orbital_or_peel_index = 1, orbital_count
            bra_occupations(orbital_or_peel_index) = csf_occupation(orbital_or_peel_index,bra_csf_index)
            ket_occupations(orbital_or_peel_index) = csf_occupation(orbital_or_peel_index,ket_csf_index)
            do quantum_number_selector = 1, 3
                bra_subshell_quantum_codes(quantum_number_selector, &
                    & orbital_or_peel_index) = csf_subshell_quantum_code(quantum_number_selector, &
                    & orbital_or_peel_index,bra_csf_index)
                ket_subshell_quantum_codes(quantum_number_selector, &
                    & orbital_or_peel_index) = csf_subshell_quantum_code(quantum_number_selector, &
                    & orbital_or_peel_index,ket_csf_index)
            end do
        end do
        ! For this CSF pair, the core consists of shells full in both states. Shells empty in both
        ! states are skipped, and all other shells enter the shared peel ordering.
        peel_orbital_count = 0
        core_orbital_count = 0
        do orbital_or_peel_index = 1, orbital_count
            if (classify_subshell_occupation(orbital_or_peel_index, &
                & bra_csf_index)==(-1) .and. classify_subshell_occupation(orbital_or_peel_index,ket_csf_index)==(-1)) cycle
            if (classify_subshell_occupation(orbital_or_peel_index, &
                & bra_csf_index)==1 .and. classify_subshell_occupation(orbital_or_peel_index,ket_csf_index)==1) then
                core_orbital_count = core_orbital_count + 1
                core_orbital_indices(core_orbital_count) = orbital_or_peel_index
            else
                peel_orbital_count = peel_orbital_count + 1
                peel_orbital_indices(peel_orbital_count) = orbital_or_peel_index
            end if
        end do
        do orbital_or_peel_index = 2, peel_orbital_count
            peel_orbital = peel_orbital_indices(orbital_or_peel_index)
            bra_intermediate_couplings(orbital_or_peel_index-1) = csf_cumulative_coupling(peel_orbital,bra_csf_index)
            ket_intermediate_couplings(orbital_or_peel_index-1) = csf_cumulative_coupling(peel_orbital,ket_csf_index)
        end do
    end subroutine prepare_csf_pair_quantum_numbers
end module csf_pair_quantum_numbers
