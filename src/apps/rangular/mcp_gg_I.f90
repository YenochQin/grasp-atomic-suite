      module mcp_i
      implicit none (type, external)
      interface
      subroutine mcp (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, output_file_prefix, temporary_file_prefix)
      implicit none (type, external)
      integer  :: symmetry_block_index
      integer  :: mpi_rank
      integer  :: mpi_process_count
      logical, intent(in) :: is_restart
      character(len=*), intent(in) :: output_file_prefix, temporary_file_prefix
      end subroutine mcp
      end interface
      end module mcp_i
