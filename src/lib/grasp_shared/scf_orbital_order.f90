module scf_orbital_order
    ! Module dependencies
    use grasp_parameters,    only: max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: orbital_update_order(max_orbitals)
end module scf_orbital_order
