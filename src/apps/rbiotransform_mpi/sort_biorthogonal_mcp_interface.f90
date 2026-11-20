module sort_biorthogonal_mcp_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine sort_biorthogonal_mcp (nfile, number, kstart, name, kamax)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nfile
            integer, intent(in) :: number
            integer :: kstart
            character (len = 24) :: name
            integer, intent(in) :: kamax
        end subroutine sort_biorthogonal_mcp
    end interface
end module sort_biorthogonal_mcp_interface
