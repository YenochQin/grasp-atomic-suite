module compute_mcp_angular_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine compute_mcp_angular_coefficients (symmetry_block_index, is_restart, mpi_rank, mpi_process_count, &
            & output_file_prefix, temporary_file_prefix)
            implicit none (type, external)

            ! Arguments
            integer  :: symmetry_block_index
            integer  :: mpi_rank
            integer  :: mpi_process_count
            logical, intent(in) :: is_restart
            character(len=*), intent(in) :: output_file_prefix, temporary_file_prefix
        end subroutine compute_mcp_angular_coefficients
    end interface
end module compute_mcp_angular_coefficients_interface
