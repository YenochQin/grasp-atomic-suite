module connect_mcp_files_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !cjb  myid, nprocs = NOT args
        subroutine connect_mcp_files (core_subshell_count, block_labels, filehead)
            implicit none (type, external)

            ! Arguments
            !     INTEGER, INTENT(IN) :: MYID
            !     INTEGER, INTENT(IN) :: NPROCS
            integer, intent(in) :: core_subshell_count
            character (len = 8), intent(in) :: block_labels(*)
            character (len = *), intent(in) :: filehead
            !...This routine performs I/O.
        end subroutine connect_mcp_files
    end interface
end module connect_mcp_files_interface
