module open_temporary_mcp_channel_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_temporary_mcp_channel (nb, kmax, filehead)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nb
            integer, intent(in) :: kmax
            character (len = *), intent(in) :: filehead
        end subroutine open_temporary_mcp_channel
    end interface
end module open_temporary_mcp_channel_interface
