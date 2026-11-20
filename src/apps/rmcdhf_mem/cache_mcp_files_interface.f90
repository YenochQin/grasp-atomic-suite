module cache_mcp_files_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !   Written by  G. Gaigalas               Vilnius, September 2021
        subroutine cache_mcp_files(nblock,nfile)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: nblock
            integer, intent(in) :: nfile
        end subroutine cache_mcp_files
    end interface
end module cache_mcp_files_interface
