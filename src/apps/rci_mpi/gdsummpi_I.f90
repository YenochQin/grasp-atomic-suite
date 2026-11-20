      module genintrk_i
      implicit none (type, external)
      interface
      subroutine genintrk (mpi_rank, mpi_process_count, n, max_orbital_two_j)
      implicit none (type, external)
      integer, intent(in) :: mpi_rank, mpi_process_count
      integer, intent(out) :: n, max_orbital_two_j
      end subroutine genintrk
      end interface
      end module genintrk_i
