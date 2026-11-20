module prepare_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine prepare_biorthogonal_mcp (name, ic, ntestg, inpci)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            integer :: ic
            integer :: ntestg
            integer :: inpci
        end subroutine prepare_biorthogonal_mcp
    end interface
end module prepare_biorthogonal_mcp_interface
