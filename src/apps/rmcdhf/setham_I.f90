      module setham_i
      implicit none (type, external)
      interface
      subroutine setham (symmetry_block_index, mpi_rank, mpi_process_count)
      implicit none (type, external)
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      end subroutine setham
      end interface
      end module setham_i
