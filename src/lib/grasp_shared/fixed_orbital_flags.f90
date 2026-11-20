module fixed_orbital_flags
    ! Module dependencies
    use grasp_parameters,    only: max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: fixed_orbital_count
    logical :: is_orbital_fixed(max_orbitals)
end module fixed_orbital_flags
