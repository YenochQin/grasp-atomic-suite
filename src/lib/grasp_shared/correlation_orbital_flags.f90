module correlation_orbital_flags
    ! Module dependencies
    use grasp_parameters,    only: max_orbitals
    implicit none (type, external)

    ! Module state
    logical :: is_correlation_orbital(max_orbitals)
end module correlation_orbital_flags
