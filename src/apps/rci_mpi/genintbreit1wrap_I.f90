      module genintbreit1wrap_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine genintbreit1wrap (mpi_rank, mpi_process_count, max_orbital_two_j)
      implicit none (type, external)
      integer, intent(in) :: mpi_rank, mpi_process_count
      integer             :: max_orbital_two_j
      end subroutine genintbreit1wrap
      end interface
      end module genintbreit1wrap_i
