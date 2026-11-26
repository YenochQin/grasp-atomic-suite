! Build the nuclear polarization potential on rank 0 and replicate it across the communicator.
module mpi_nuclear_polarization
    ! Module dependencies
    use mpi_f08, only: mpi_comm, mpi_comm_rank, mpi_bcast, mpi_allreduce, mpi_abort, &
        mpi_integer, mpi_double_precision, mpi_logical, mpi_character, mpi_min
    use grasp_parameters, only: max_radial_points
    use nuclear_polarization_state, only: npol_capture, npol_config, npol_message_len, npol_publish, npol_reset, &
        coulomb_potential_baseline
    use nuclear_polarization_input, only: npol_read_config, npol_validate_config
    use nuclear_polarization_potential, only: npol_check_grid, npol_build_grid, npol_compose
    use nuclear_polarization_initialization, only: npol_write_summary, npol_summary_path
    implicit none (type, external)

    ! Module interface
    private
    public :: npol_init_mpi, npol_abort_mpi
contains
    subroutine npol_init_mpi(config_path, nuclear_charge, grid_point_count, radial_grid, negative_radius_times_potential, &
        communicator, status, message)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64

        ! Arguments
        character(*), intent(in) :: config_path
        real(real64), intent(in) :: nuclear_charge, radial_grid(:)
        integer, intent(in) :: grid_point_count
        type(mpi_comm), intent(in) :: communicator
        real(real64), intent(inout) :: negative_radius_times_potential(:)
        integer, intent(out) :: status
        character(*), intent(out) :: message

        ! Local declarations
        type(npol_config) :: config
        real(real64) :: dipole_potential(max_radial_points), quadrupole_potential(max_radial_points)
        real(real64) :: polarization_potential(max_radial_points), composed_potential(max_radial_points)
        integer :: rank, mpi_return_code, root_grid_point_count
        ! Use a fixed message length for status broadcasts.
        character(npol_message_len) :: status_message
        call mpi_comm_rank(communicator, rank, mpi_return_code)
        call npol_reset()
        call npol_check_grid(grid_point_count, radial_grid, status, status_message)
        if (size(negative_radius_times_potential) < grid_point_count) then
            status = 1
            status_message = "npol: Coulomb array shorter than N"
        end if
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        root_grid_point_count = grid_point_count
        call mpi_bcast(root_grid_point_count, 1, mpi_integer, 0, communicator, mpi_return_code)
        if (root_grid_point_count /= grid_point_count) then
            status = 1
            status_message = "npol: inconsistent active grid length across ranks"
        end if
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        call npol_capture(grid_point_count, negative_radius_times_potential(:grid_point_count))

        if (rank == 0) then
            call npol_read_config(config_path, config, status, status_message)
            if (status == 0) call npol_validate_config(config, nuclear_charge, status, status_message)
        end if
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        call broadcast_config(config, communicator)
        ! Validate the shared configuration against each rank's nuclear charge.
        call npol_validate_config(config, nuclear_charge, status, status_message)
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        if (rank == 0) call npol_build_grid(config, grid_point_count, radial_grid, dipole_potential, quadrupole_potential, &
            polarization_potential, status, status_message)
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        call mpi_bcast(dipole_potential, grid_point_count, mpi_double_precision, 0, communicator, mpi_return_code)
        call mpi_bcast(quadrupole_potential, grid_point_count, mpi_double_precision, 0, communicator, mpi_return_code)
        call mpi_bcast(polarization_potential, grid_point_count, mpi_double_precision, 0, communicator, mpi_return_code)
        call npol_compose(grid_point_count, radial_grid, coulomb_potential_baseline, polarization_potential, &
            composed_potential, status, status_message)
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        if (rank == 0) call npol_write_summary(npol_summary_path(config_path), config, grid_point_count, status, status_message)
        call sync_status(communicator, rank, status, status_message)
        message = status_message
        if (status /= 0) return
        ! Publish the full -rV grid after all ranks pass the composition and summary checks.
        ! YPOT apportions the external potential during the later COFPOT reduction.
        negative_radius_times_potential(:grid_point_count) = composed_potential(:grid_point_count)
        call npol_publish(config, grid_point_count, dipole_potential(:grid_point_count), &
            quadrupole_potential(:grid_point_count), polarization_potential(:grid_point_count))
    end subroutine npol_init_mpi

    subroutine sync_status(communicator, rank, status, message)
        ! Arguments
        type(mpi_comm), intent(in) :: communicator
        integer, intent(in) :: rank
        integer, intent(inout) :: status
        character(npol_message_len), intent(inout) :: message

        ! Local declarations
        integer :: failed_rank, first_failed_rank, mpi_return_code
        ! The minimum failing rank supplies the error message for the communicator.
        failed_rank = huge(0)
        if (status /= 0) failed_rank = rank
        call mpi_allreduce(failed_rank, first_failed_rank, 1, mpi_integer, mpi_min, communicator, mpi_return_code)
        status = 0
        if (first_failed_rank == huge(0)) then
            message = ""
            return
        end if
        status = 1
        call mpi_bcast(message, len(message), mpi_character, first_failed_rank, communicator, mpi_return_code)
    end subroutine sync_status

    subroutine broadcast_config(config, communicator)
        ! Arguments
        type(npol_config), intent(inout) :: config
        type(mpi_comm), intent(in) :: communicator

        ! Local declarations
        integer :: mpi_return_code
        ! Broadcast fields with their MPI types so the transfer is independent of derived-type padding.
        call mpi_bcast(config%schema_version, 1, mpi_integer, 0, communicator, mpi_return_code)
        call mpi_bcast(config%isotope_z, 1, mpi_integer, 0, communicator, mpi_return_code)
        call mpi_bcast(config%isotope_a, 1, mpi_integer, 0, communicator, mpi_return_code)
        call mpi_bcast(config%enabled, 1, mpi_logical, 0, communicator, mpi_return_code)
        call mpi_bcast(config%use_e1, 1, mpi_logical, 0, communicator, mpi_return_code)
        call mpi_bcast(config%use_e2, 1, mpi_logical, 0, communicator, mpi_return_code)
        call mpi_bcast(config%k1, 1, mpi_double_precision, 0, communicator, mpi_return_code)
        call mpi_bcast(config%b1, 1, mpi_double_precision, 0, communicator, mpi_return_code)
        call mpi_bcast(config%k2, 1, mpi_double_precision, 0, communicator, mpi_return_code)
        call mpi_bcast(config%b2, 1, mpi_double_precision, 0, communicator, mpi_return_code)
        call mpi_bcast(config%units, len(config%units), mpi_character, 0, communicator, mpi_return_code)
        call mpi_bcast(config%model_id, len(config%model_id), mpi_character, 0, communicator, mpi_return_code)
        call mpi_bcast(config%parameter_source, len(config%parameter_source), mpi_character, 0, communicator, mpi_return_code)
    end subroutine broadcast_config

    subroutine npol_abort_mpi(communicator, message)
        ! Arguments
        type(mpi_comm), intent(in) :: communicator
        character(*), intent(in) :: message

        ! Local declarations
        integer :: rank, mpi_return_code
        call mpi_comm_rank(communicator, rank, mpi_return_code)
        if (rank == 0) write(*,"(A)") trim(message)
        call mpi_abort(communicator, 1, mpi_return_code)
    end subroutine npol_abort_mpi
end module mpi_nuclear_polarization
