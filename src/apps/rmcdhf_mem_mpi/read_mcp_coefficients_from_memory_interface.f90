module read_mcp_coefficients_from_memory_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !   Written by  G. Gaigalas               Vilnius, September 2021
        subroutine read_mcp_coefficients_from_memory(jblock)
            implicit none (type, external)

            ! Arguments
            integer, intent(in)       :: jblock
        end subroutine read_mcp_coefficients_from_memory
    end interface
end module read_mcp_coefficients_from_memory_interface
