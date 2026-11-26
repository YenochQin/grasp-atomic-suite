! Expose mpi_f08 and the suite's shared process data.
module mpi_state
    ! Module dependencies
    use mpi_f08
    use mpi_runtime, only: initialize_mpi_runtime, create_contiguous_mpi_byte_datatype
    implicit none (type, external)

    ! Module state
    character(len=mpi_max_processor_name) :: processor_name
    integer :: mpi_rank, mpi_process_count, processor_name_length
    integer :: mpi_error_code

end module mpi_state
