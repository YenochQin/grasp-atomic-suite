module biorthogonal_subshell_mapping
    ! Module dependencies
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module constants
    integer, parameter :: nlmax = 20

    ! Module state
    integer :: nshlp(nlmax, nlmax)
    integer :: nshlpp(nlmax, max_orbitals)
end module biorthogonal_subshell_mapping
