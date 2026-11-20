module angular_coupling_workspace
    implicit none (type, external)

    ! Module constants
    integer, parameter :: angular_value_capacity = 60
    integer, parameter :: coupling_triad_capacity = 12

    ! Module state
    integer :: coupling_two_j_plus_one(angular_value_capacity)
    integer :: bra_coupling_triads(coupling_triad_capacity, 3), ket_coupling_triads(coupling_triad_capacity, 3)
    integer :: coupling_angular_value_count, coupling_triad_count
    logical :: free_coupling_nodes(angular_value_capacity)
end module angular_coupling_workspace
