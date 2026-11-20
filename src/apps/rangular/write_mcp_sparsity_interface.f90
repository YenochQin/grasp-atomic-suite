module write_mcp_sparsity_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_mcp_sparsity (write_coefficient_callback, nonzero_element_count, print_details, symmetry_block_index, &
            & mpi_rank, &
            & mpi_process_count, output_file_prefix)
            ! Procedure interfaces
            use write_mcp_sparsity_summary_interface, only: output_callback => write_mcp_sparsity_summary
            implicit none (type, external)

            ! Arguments
            procedure(output_callback) :: write_coefficient_callback
            integer, intent(in) :: nonzero_element_count
            logical :: print_details
            integer, intent(in) :: symmetry_block_index
            integer, intent(in) :: mpi_rank
            integer, intent(in) :: mpi_process_count
            character (len = *), intent(in) :: output_file_prefix
        end subroutine write_mcp_sparsity
    end interface
end module write_mcp_sparsity_interface
