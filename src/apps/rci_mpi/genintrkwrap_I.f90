      module genintrkwrap_i
      implicit none (type, external)
      interface
      subroutine genintrkwrap (mpi_rank, mpi_process_count, max_orbital_two_j)
      implicit none (type, external)
      integer :: mpi_rank
      integer :: mpi_process_count
      integer :: max_orbital_two_j
      end subroutine genintrkwrap
      end interface
      end module genintrkwrap_i
