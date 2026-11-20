module write_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_biorthogonal_mcp (name, ik, ntestg, inpci)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            integer, intent(in) :: ik
            integer, intent(in) :: ntestg
            integer :: inpci
        end subroutine write_biorthogonal_mcp
    end interface
end module write_biorthogonal_mcp_interface
