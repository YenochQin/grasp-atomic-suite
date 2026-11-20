      module genintbreit2_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine genintbreit2 (mpi_rank, mpi_process_count, nb, max_orbital_two_j)
      implicit none (type, external)
      integer, intent(in) :: mpi_rank, mpi_process_count
      integer             :: nb, max_orbital_two_j
      end subroutine genintbreit2
      end interface
      end module genintbreit2_i
