module angular_recoupling_backup
    ! Module dependencies
    use numeric_constants, only: eleven, angular_tolerance, four, half, one, seven, ten, tenth, three, two, zero
    use grasp_parameters,   only:  max_orbitals
    implicit none (type, external)

    ! Module state
    integer :: saved_peel_orbitals(max_orbitals), saved_bra_j_weights(max_orbitals), saved_ket_j_weights(max_orbitals)
end module angular_recoupling_backup
