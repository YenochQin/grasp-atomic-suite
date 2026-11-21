! Expose the MPI Fortran interfaces and the suite's shared rank, processor and status data.
module mpi_state
    ! Module dependencies
    use mpi
    use mpi_runtime, only: initialize_mpi_runtime, create_contiguous_mpi_byte_datatype
    implicit none (type, external)

    ! Module state
    character(len=mpi_max_processor_name) :: processor_name
    integer :: mpi_rank, mpi_process_count, processor_name_length
    integer :: mpi_error_code
    integer :: mpi_status_buffer(mpi_status_size)

end module mpi_state
