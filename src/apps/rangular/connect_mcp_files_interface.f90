module connect_mcp_files_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine connect_mcp_files (mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: mpi_rank
            integer, intent(in) :: mpi_process_count
            integer, intent(in) :: core_subshell_count
            character (len = 8), intent(in) :: block_labels(*)
            character (len = *), intent(in) :: filehead
            !...This routine performs I/O.
        end subroutine connect_mcp_files
    end interface
end module connect_mcp_files_interface
