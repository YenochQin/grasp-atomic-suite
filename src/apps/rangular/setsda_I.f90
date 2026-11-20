      module setsda_i
      implicit none (type, external)
      interface
      subroutine setsda (write_coefficient_callback, nonzero_element_count, print_details, symmetry_block_index, mpi_rank, &
         & mpi_process_count, output_file_prefix)
      use outsda_i, only: output_callback => outsda
      implicit none (type, external)
      procedure(output_callback) :: write_coefficient_callback
      integer, intent(in) :: nonzero_element_count
      logical :: print_details
      integer, intent(in) :: symmetry_block_index
      integer, intent(in) :: mpi_rank
      integer, intent(in) :: mpi_process_count
      character (len = *), intent(in) :: output_file_prefix
      end subroutine setsda
      end interface
      end module setsda_i
