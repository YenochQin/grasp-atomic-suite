! Select ASF roots and open the indexed CI mixing stream.
module open_mixing_output_mpi_interface
    implicit none (type, external)
    private
    public :: open_mixing_output_mpi
contains
    subroutine open_mixing_output_mpi(name, block_labels)
        use calculation_state_files, only: ci_mixing_suffix
        use calculation_file_units, only: mixing_output_unit
        use read_ci_mixing_coefficients_mpi_interface, only: read_ci_mixing_coefficients_mpi
        use mixing_file_format, only: open_mixing_stream
        use mpi_state, only: mpi_rank
        use mpi_io_errors, only: fail_input_io_mpi
        character(*), intent(in) :: name
        character(8), intent(in) :: block_labels(*)
        integer :: status

        call read_ci_mixing_coefficients_mpi(block_labels)
        if (mpi_rank /= 0) return
        call open_mixing_stream(trim(name)//ci_mixing_suffix, mixing_output_unit, status)
        if (status /= 0) call fail_input_io_mpi(status, mixing_output_unit, "Cannot create CI mixing stream")
    end subroutine
end module open_mixing_output_mpi_interface
