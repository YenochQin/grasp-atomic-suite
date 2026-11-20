      module genintbreit2wrap_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine genintbreit2wrap (mpi_rank, mpi_process_count, max_orbital_two_j)
      implicit none (type, external)
      integer, intent(in) :: mpi_rank, mpi_process_count
      integer             :: max_orbital_two_j
      end subroutine genintbreit2wrap
      end interface
      end module genintbreit2wrap_i
