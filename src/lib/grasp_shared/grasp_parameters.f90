! Re-exports the central capacities and orbital packing radix.
module grasp_parameters

    ! Module dependencies
    use compile_time_capacities, only: orbital_index_radix, max_radial_points, radial_workspace_capacity, &
        max_orbitals, max_orbitals_minus_one, max_orbitals_minus_two, csf_total_j_slot, csf_parity_slot, csf_angular_code_count, &
        csf_orbital_letters, maximum_orbital_two_j
    implicit none (type, external)
end module grasp_parameters
