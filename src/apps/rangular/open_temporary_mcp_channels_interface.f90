module open_temporary_mcp_channels_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine open_temporary_mcp_channels (nb, k, filehead)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nb
            integer, intent(in) :: k
            character (len = *), intent(in) :: filehead
        end subroutine open_temporary_mcp_channels
    end interface
end module open_temporary_mcp_channels_interface
