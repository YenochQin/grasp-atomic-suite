      module setmcp2_i
      implicit none (type, external)
      interface
      subroutine setmcp2 (mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
      implicit none (type, external)
      integer, intent(in) :: mpi_rank
      integer :: mpi_process_count
      integer :: core_subshell_count
      character (len = 8) :: block_labels(*)
      character (len = *) :: filehead
      end subroutine setmcp2
      end interface
      end module setmcp2_i
