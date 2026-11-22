! Supply the rank and process count used by the serial variants of the MCDHF kernels.
module serial_mpi_context
    implicit none (type, external)

    ! Module state
    integer :: mpi_rank, mpi_process_count
end module serial_mpi_context
