! Report input errors with an MPI abort callback for an active MPI session.
module mpi_io_errors
    ! Module dependencies
    use mpi_f08, only: mpi_comm_world, mpi_initialized, mpi_finalized, mpi_abort
    use input_error_handling, only: fail_input_io
    implicit none (type, external)

    ! Module interface
    private
    public :: fail_input_io_mpi, abort_input_mpi

contains

    subroutine fail_input_io_mpi(io_status, io_unit, context)
        ! Arguments
        integer, intent(in) :: io_status, io_unit
        character(*), intent(in) :: context
        call fail_input_io(io_status, io_unit, context, abort_input_mpi)
    end subroutine fail_input_io_mpi

    subroutine abort_input_mpi()
        ! Local declarations
        integer :: mpi_return_code
        logical :: mpi_is_initialized, mpi_is_finalized
        call mpi_initialized(mpi_is_initialized, mpi_return_code)
        if (.not. mpi_is_initialized) return
        call mpi_finalized(mpi_is_finalized, mpi_return_code)
        if (.not. mpi_is_finalized) call mpi_abort(mpi_comm_world, 1, mpi_return_code)
    end subroutine abort_input_mpi

end module mpi_io_errors
