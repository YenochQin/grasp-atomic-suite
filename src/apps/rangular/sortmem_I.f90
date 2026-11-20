      module sortmem_i
      implicit none (type, external)
      interface
      subroutine sortmem (mcp_channel_index,sort_capacity_per_channel,coefficient_count,distinct_integral_count,print_details, &
         & symmetry_block_index,output_file_prefix)
      implicit none (type, external)
      integer, intent(in)  :: mcp_channel_index
      integer, intent(in)  :: coefficient_count
      integer, intent(out) :: distinct_integral_count
      logical, intent(in)  :: print_details
      integer, intent(in)  :: symmetry_block_index,sort_capacity_per_channel
      character (len = *), intent(in) :: output_file_prefix
      end subroutine sortmem
      end interface
      end module sortmem_i
