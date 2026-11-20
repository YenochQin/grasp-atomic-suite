module initial_orbital_source
    ! Module dependencies
    use grasp_parameters,   only: max_orbitals
    implicit none (type, external)

    ! Module state
    character :: orbital_source_description(max_orbitals)*256
end module initial_orbital_source
