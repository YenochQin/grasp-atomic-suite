      module mpi_c
      use, intrinsic :: iso_fortran_env, only: double => real64
! Use the MPI implementation's explicit Fortran interfaces. Integer handles
! remain the public convention of this legacy module; choice buffers are
! described by the MPI module rather than unchecked external declarations.
      use mpi
      use mpi_runtime_i, only: startmpi, mpix_bytes
      implicit none (type, external)

      character(len=mpi_max_processor_name) :: processor_name
      integer :: mpi_rank, mpi_process_count, processor_name_length
      integer :: mpi_error_code, ierrtotal
      integer :: mpi_status_buffer(mpi_status_size)



      end module mpi_c
