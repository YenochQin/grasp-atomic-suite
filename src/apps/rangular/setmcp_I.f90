      module setmcp_i
      implicit none (type, external)
      interface
      subroutine setmcp (mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
      implicit none (type, external)
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      integer, intent(in) :: core_subshell_count
      character (len = 8), intent(in) :: block_labels(*)
      character (len = *), intent(in) :: filehead
!...This routine performs I/O.
      end subroutine setmcp
      end interface
      end module setmcp_i
