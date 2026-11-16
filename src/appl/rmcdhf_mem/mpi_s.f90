      module mpi_s
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer :: mpi_rank, mpi_process_count
      end module mpi_s
