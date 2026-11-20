subroutine build_orbital_potential(extended_optimal_level, orbital_index, radial_point_count)
    ! Module dependencies
    !-----------------------------------------------
    use scf_radial_potentials, only: max_radial_points, xp, xq, yp
    use serial_mpi_context, only: mpi_process_count

    ! Procedure interfaces
    use assemble_orbital_potential_coefficients_interface, only: assemble_orbital_potential_coefficients
    use build_direct_radial_potential_interface, only: build_direct_radial_potential
    use build_exchange_radial_potential_interface, only: build_exchange_radial_potential
    use apply_lagrange_constraints_interface, only: apply_lagrange_constraints
    use apply_off_diagonal_integrals_interface, only: apply_off_diagonal_integrals
    implicit none (type, external)

    ! Arguments
    integer  :: orbital_index
    integer  :: radial_point_count
    logical  :: extended_optimal_level
    !-----------------------------------------------------------------------

    call assemble_orbital_potential_coefficients (extended_optimal_level, orbital_index)
    call build_direct_radial_potential (orbital_index)
    call build_exchange_radial_potential (orbital_index)
    call apply_lagrange_constraints (orbital_index, mpi_process_count)
    call apply_off_diagonal_integrals

    return
end subroutine build_orbital_potential
