module sort_mcp_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine sort_mcp_coefficients (mcp_channel_index, coefficient_count, distinct_integral_count, print_details, &
            & symmetry_block_index, &
            & output_file_prefix)
            implicit none (type, external)

            ! Arguments
            integer, intent(in)  :: mcp_channel_index
            integer, intent(in)  :: coefficient_count
            integer, intent(out) :: distinct_integral_count
            logical, intent(in)  :: print_details
            integer, intent(in)  :: symmetry_block_index
            character (len = *), intent(in) :: output_file_prefix
        end subroutine sort_mcp_coefficients
    end interface
end module sort_mcp_coefficients_interface
