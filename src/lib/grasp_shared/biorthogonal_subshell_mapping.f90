module biorthogonal_subshell_mapping
    ! Module dependencies
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module constants
    integer, parameter :: active_kappa_group_capacity = 20

    ! Module state
    integer :: active_orbitals_by_group(active_kappa_group_capacity, active_kappa_group_capacity)
    integer :: active_orbital_positions_by_group(active_kappa_group_capacity, max_orbitals)
end module biorthogonal_subshell_mapping
