module read_mcp_coefficients_into_memory_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !   Written by  G. Gaigalas               Vilnius, September 2021
        subroutine read_mcp_coefficients_into_memory(nfile,total_contribution_count,ncontr)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: int64
            implicit none (type, external)

            ! Arguments
            integer(int64), intent(in) :: total_contribution_count
            integer, intent(in)       :: nfile
            integer, intent(out)      :: ncontr
        end subroutine read_mcp_coefficients_into_memory
    end interface
end module read_mcp_coefficients_into_memory_interface
