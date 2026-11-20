module read_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine read_biorthogonal_mcp (name, startup_directory, ik, ntestg, inpci)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
            character(len=128), intent(in) :: startup_directory
            integer, intent(in) :: ik
            integer, intent(in) :: ntestg
            integer, intent(in) :: inpci
        end subroutine read_biorthogonal_mcp
    end interface
end module read_biorthogonal_mcp_interface
