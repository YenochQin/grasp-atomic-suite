module read_mcp_structure_into_memory_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !   Written by  G. Gaigalas               Vilnius, September 2021
        subroutine read_mcp_structure_into_memory(nfile,total_contribution_count,ncontr,iclmndum)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: int64
            implicit none (type, external)

            ! Arguments
            integer(int64), intent(in) :: total_contribution_count
            integer, intent(in)       :: nfile
            integer, intent(out)      :: ncontr
            integer, intent(out)      :: iclmndum
        end subroutine read_mcp_structure_into_memory
    end interface
end module read_mcp_structure_into_memory_interface
