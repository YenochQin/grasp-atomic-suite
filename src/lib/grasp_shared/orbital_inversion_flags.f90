module orbital_inversion_flags
    ! Module dependencies
    use grasp_parameters,    only: max_orbitals
    implicit none (type, external)

    ! Module state
    logical :: noinvt(max_orbitals)
end module orbital_inversion_flags
