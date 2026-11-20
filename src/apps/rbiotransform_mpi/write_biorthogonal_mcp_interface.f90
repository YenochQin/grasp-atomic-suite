module write_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_biorthogonal_mcp (name, startup_directory,ik, ntestg, inpci)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            character (len = 128) :: startup_directory
            integer, intent(in) :: ik
            integer, intent(in) :: ntestg
            integer :: inpci
        end subroutine write_biorthogonal_mcp
    end interface
end module write_biorthogonal_mcp_interface
