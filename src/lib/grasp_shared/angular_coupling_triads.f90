module angular_coupling_triads
    implicit none (type, external)

    ! Module constants
    integer, parameter :: saved_coupling_triad_capacity = 12

    ! Module state
    integer :: saved_bra_coupling_triads(saved_coupling_triad_capacity, 3), &
        saved_ket_coupling_triads(saved_coupling_triad_capacity, 3)
end module angular_coupling_triads
