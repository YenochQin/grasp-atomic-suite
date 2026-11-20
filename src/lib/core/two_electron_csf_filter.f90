! Reject CSF pairs whose occupation differences exceed the two-electron screening bounds.
module two_electron_csf_filter
    implicit none (type, external)

    ! Module interface
    private
    public :: ichkq2
contains
    integer function ichkq2 (bra_csf_index, ket_csf_index)
        ! Module dependencies
        use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
        use orbital_csf_basis, only: orbital_count
        use csf_quantum_queries, only: iq
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index
        integer, intent(in) :: ket_csf_index

        ! Local declarations
        integer :: differing_subshell_count, orbital_index, bra_occupation, ket_occupation
        ichkq2 = 0
        differing_subshell_count = 0
        do orbital_index = 1, orbital_count
            bra_occupation = iq(orbital_index,bra_csf_index)
            ket_occupation = iq(orbital_index,ket_csf_index)
            if (bra_occupation == ket_occupation) cycle
            differing_subshell_count = differing_subshell_count + 1
            if (differing_subshell_count > 4) return
            if (iabs(bra_occupation - ket_occupation) <= 2) cycle
            return
        end do
        ichkq2 = 1
        return
    end function ichkq2
end module two_electron_csf_filter
