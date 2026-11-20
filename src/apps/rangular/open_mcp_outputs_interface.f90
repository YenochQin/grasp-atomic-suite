module open_mcp_outputs_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_mcp_outputs (mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: mpi_rank
            integer :: mpi_process_count
            integer :: core_subshell_count
            character (len = 8) :: block_labels(*)
            character (len = *) :: filehead
        end subroutine open_mcp_outputs
    end interface
end module open_mcp_outputs_interface
