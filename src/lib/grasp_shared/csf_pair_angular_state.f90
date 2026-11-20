! Holds occupations, coupling data, and the peel/core orbital lists for the current CSF pair.
module csf_pair_angular_state
    ! Module dependencies
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: nq1(max_orbitals), nq2(max_orbitals)
    ! setqna builds these chains in peel order, with one coupling node per shell after the first.
    integer :: jjc1(max_orbitals), jjc2(max_orbitals)
    ! These pair-local arrays put the quantum-number selector first for each orbital.
    integer :: jjq1(3, max_orbitals), jjq2(3, max_orbitals)
    integer :: jlist(max_orbitals), klist(max_orbitals)
    integer :: npeel, ncore
end module csf_pair_angular_state
