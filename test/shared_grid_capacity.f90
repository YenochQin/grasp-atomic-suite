! Use an independent module so the expected values come directly from the central definitions.
module shared_grid_capacity_reference
    use cal_constants, only: max_radial_points, radial_workspace_capacity, orbital_index_radix, max_orbitals, &
        max_orbitals_minus_one, max_orbitals_minus_two
    implicit none (type, external)
    private
    integer, parameter, public :: reference_capacities(*) = &
        [max_radial_points, radial_workspace_capacity, orbital_index_radix, max_orbitals, &
         max_orbitals_minus_one, max_orbitals_minus_two]
end module shared_grid_capacity_reference

program shared_grid_capacity
    use shared_grid_capacity_reference, only: reference_capacities
    use grasp_parameters, only: max_radial_points, radial_workspace_capacity, orbital_index_radix, max_orbitals, &
        max_orbitals_minus_one, max_orbitals_minus_two, csf_total_j_slot, csf_parity_slot, csf_angular_code_count
    use radial_grid_state, only: radial_grid, radial_grid_derivative, radial_derivative_over_radius
    use scf_radial_potentials, only: large_exchange_source, small_exchange_source, orbital_y_potential
    use radial_quadrature_workspace, only: radial_workspace_a, radial_workspace_b
    implicit none (type, external)

    if (any([max_radial_points, radial_workspace_capacity, orbital_index_radix, max_orbitals, &
             max_orbitals_minus_one, max_orbitals_minus_two] /= reference_capacities)) &
        error stop "grasp_parameters does not re-export the configured limits"
    if (radial_workspace_capacity /= max_radial_points + 10) error stop "inconsistent radial capacity parameters"
    if (csf_total_j_slot /= max_orbitals+1 .or. csf_parity_slot /= max_orbitals+2 .or. &
        csf_angular_code_count /= csf_parity_slot) error stop "inconsistent CSF angular code capacity"
    if (max_orbitals_minus_one /= max_orbitals - 1 .or. max_orbitals_minus_two /= max_orbitals - 2) &
        error stop "inconsistent orbital capacities"
    if (any([size(radial_grid), size(radial_grid_derivative), size(radial_derivative_over_radius)] &
            /= radial_workspace_capacity)) error stop "radial_grid_state was compiled with a different radial capacity"
    if (any([size(large_exchange_source), size(small_exchange_source), size(orbital_y_potential)] /= max_radial_points)) &
        error stop "scf_radial_potentials was compiled with a different radial capacity"
    if (size(radial_workspace_a) /= radial_workspace_capacity .or. size(radial_workspace_b) /= radial_workspace_capacity) &
        error stop "radial_quadrature_workspace was compiled with a different radial capacity"
    print *, "Shared radial capacity:", max_radial_points, radial_workspace_capacity
end program shared_grid_capacity
