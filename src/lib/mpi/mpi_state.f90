      ! Expose the MPI Fortran interfaces and the suite's shared rank, processor and status data.
      module mpi_state
      use mpi
      use mpi_runtime, only: startmpi, mpix_bytes
      implicit none (type, external)

      character(len=mpi_max_processor_name) :: processor_name
      integer :: mpi_rank, mpi_process_count, processor_name_length
      integer :: mpi_error_code, ierrtotal
      integer :: mpi_status_buffer(mpi_status_size)

      end module mpi_state
