module open_mcp_outputs_mpi_interface
    ! Module dependencies
    use csf_descriptor_reader, only: descriptor_source
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !cjb  myid, nprocs = NOT args
        subroutine open_mcp_outputs_mpi (core_subshell_count, block_labels, filehead, descriptors)
            ! Module dependencies
            import descriptor_source
            implicit none (type, external)

            ! Arguments
            !     INTEGER , INTENT(IN) :: MYID
            !     INTEGER , INTENT(IN) :: NPROCS
            integer , intent(in) :: core_subshell_count
            character (len = 8), intent(in) :: block_labels(*)
            character (len = *), intent(in) :: filehead
            type(descriptor_source), intent(in) :: descriptors
            !...This routine performs I/O.
        end subroutine open_mcp_outputs_mpi
    end interface
end module open_mcp_outputs_mpi_interface
