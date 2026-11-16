module input_io_mpi
    use mpi, only: mpi_comm_world, mpi_initialized, mpi_finalized, mpi_abort
    use input_io, only: fail_input_io
    implicit none
    private
    public :: fail_input_io_mpi

contains

    subroutine fail_input_io_mpi(status, unit, context)
        integer, intent(in) :: status, unit
        character(*), intent(in) :: context
        call fail_input_io(status, unit, context, abort_input_mpi)
    end subroutine fail_input_io_mpi

    subroutine abort_input_mpi()
        integer :: status
        logical :: initialized, finalized
        call mpi_initialized(initialized, status)
        if (.not. initialized) return
        call mpi_finalized(finalized, status)
        if (.not. finalized) call mpi_abort(mpi_comm_world, 1, status)
    end subroutine abort_input_mpi

end module input_io_mpi
