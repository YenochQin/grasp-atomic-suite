! Check whether CSF occupations can be connected by a one-electron operator.
module one_electron_csf_filter
    implicit none (type, external)

    ! Module interface
    private
    public :: check_one_electron_csf_pair
contains
    integer function check_one_electron_csf_pair (bra_csf_index, ket_csf_index)
        ! Module dependencies
        use orbital_csf_basis, only: orbital_count
        use csf_quantum_queries, only: csf_occupation
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index
        integer, intent(in) :: ket_csf_index

        ! Local declarations
        integer :: differing_subshell_count, orbital_index, bra_occupation, ket_occupation

        check_one_electron_csf_pair = 0
        differing_subshell_count = 0
        do orbital_index = 1, orbital_count
            bra_occupation = csf_occupation(orbital_index,bra_csf_index)
            ket_occupation = csf_occupation(orbital_index,ket_csf_index)
            if (bra_occupation == ket_occupation) cycle
            differing_subshell_count = differing_subshell_count + 1
            if (differing_subshell_count > 2) return
            if (iabs(bra_occupation - ket_occupation) <= 1) cycle
            return
        end do
        if (differing_subshell_count==2 .or. differing_subshell_count==0) check_one_electron_csf_pair = 1
        return
    end function check_one_electron_csf_pair
end module one_electron_csf_filter
