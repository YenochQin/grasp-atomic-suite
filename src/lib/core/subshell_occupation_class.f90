! Classify a CSF subshell as empty (-1), partly occupied (0), or full (+1).
module subshell_occupation_class
    implicit none (type, external)

    ! Module interface
    private
    public :: ichop
contains
    integer function ichop (subshell_index, csf_index)
        ! Module dependencies
        use grasp_parameters, only: max_orbitals
        use orbital_csf_basis,         only: orbital_l, orbital_two_j
        use csf_quantum_queries, only: iq
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: subshell_index
        integer, intent(in) :: csf_index

        ! Local declarations
        integer :: occupation, subshell_capacity
        occupation = iq(subshell_index,csf_index)
        subshell_capacity = orbital_two_j(subshell_index) + 1
        if (occupation == 0) then
            ichop = -1
        else if (occupation == subshell_capacity) then
            ichop = 1
        else
            ichop = 0
        end if
        return
    end function ichop
end module subshell_occupation_class
