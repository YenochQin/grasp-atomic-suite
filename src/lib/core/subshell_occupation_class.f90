! Classify a CSF subshell as empty (-1), partly occupied (0), or full (+1).
module subshell_occupation_class
    implicit none (type, external)

    ! Module interface
    private
    public :: classify_subshell_occupation
contains
    integer function classify_subshell_occupation (subshell_index, csf_index)
        ! Module dependencies
        use orbital_csf_basis, only: orbital_two_j
        use csf_quantum_queries, only: csf_occupation
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: subshell_index
        integer, intent(in) :: csf_index

        ! Local declarations
        integer :: occupation, subshell_capacity
        occupation = csf_occupation(subshell_index,csf_index)
        subshell_capacity = orbital_two_j(subshell_index) + 1
        if (occupation == 0) then
            classify_subshell_occupation = -1
        else if (occupation == subshell_capacity) then
            classify_subshell_occupation = 1
        else
            classify_subshell_occupation = 0
        end if
        return
    end function classify_subshell_occupation
end module subshell_occupation_class
