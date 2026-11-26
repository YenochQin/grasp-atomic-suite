! Restore the CI calculation controls and radial state from the restart file.
module read_ci_restart_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_ci_restart
contains

    subroutine read_ci_restart
        use calculation_file_units, only: calculation_log_unit
        use hamiltonian_correction_flags, only: estimate_self_energy
        use ci_restart_units, only: restart_matrix_unit
        use standard_io_units, only: stderr_unit
        use yes_no_input, only: read_yes_no
        use rci_restart_records, only: read_rci_restart_header
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
        use mpi_io_errors, only: abort_input_mpi
        implicit none (type, external)

        call read_rci_restart_header(restart_matrix_unit, abort_input_mpi)
        ! Determine if the self-energy contribution is to be estimated
        if (mpi_rank == 0) then
            write (stderr_unit, *) "Estimate contributions from the self-energy?"
            estimate_self_energy = read_yes_no()
            if (estimate_self_energy) then
                write(calculation_log_unit,"(a)") "y            ! Estimate contributions from the self-energy?"
            else
                write(calculation_log_unit,"(a)") "n            ! Estimate contributions from the self-energy?"
            end if
        end if
        call mpi_bcast (estimate_self_energy, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
    end subroutine read_ci_restart
end module read_ci_restart_interface
