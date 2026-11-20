module read_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine read_biorthogonal_mcp (name, ik, ntestg, inpci)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            integer, intent(in) :: ik
            integer, intent(in) :: ntestg
            integer, intent(in) :: inpci
        end subroutine read_biorthogonal_mcp
    end interface
end module read_biorthogonal_mcp_interface
