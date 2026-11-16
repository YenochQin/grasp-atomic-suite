      module mcp_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine mcp (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, output_file_prefix)
      integer  :: symmetry_block_index
      integer  :: mpi_rank
      integer  :: mpi_process_count
      logical, intent(in) :: is_restart
      character(len=*), intent(in) :: output_file_prefix
      end subroutine
      end interface
      end module
