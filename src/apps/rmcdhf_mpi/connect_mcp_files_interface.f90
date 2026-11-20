module connect_mcp_files_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine connect_mcp_files (core_subshell_count, input_symmetry_block_count, block_labels, filehead)
            implicit none (type, external)

            ! Arguments
            integer :: core_subshell_count
            integer, intent(in) :: input_symmetry_block_count
            character (len = 8) :: block_labels(*)
            character (len = *), intent(in) :: filehead
            !...This routine performs I/O.
        end subroutine connect_mcp_files
    end interface
end module connect_mcp_files_interface
