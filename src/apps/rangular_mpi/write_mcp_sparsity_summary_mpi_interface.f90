module write_mcp_sparsity_summary_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_mcp_sparsity_summary_mpi (print_details, nonzero_element_count, csf_count)
            implicit none (type, external)

            ! Arguments
            logical, intent(in) :: print_details
            integer, intent(in) :: nonzero_element_count
            integer, intent(in) :: csf_count
        end subroutine write_mcp_sparsity_summary_mpi
    end interface
end module write_mcp_sparsity_summary_mpi_interface
