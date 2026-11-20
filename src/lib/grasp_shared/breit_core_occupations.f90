! Stores core occupations used by the Breit matrix-element construction.
module breit_core_occupations
    ! Module dependencies
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: icore(max_orbitals)
end module breit_core_occupations
