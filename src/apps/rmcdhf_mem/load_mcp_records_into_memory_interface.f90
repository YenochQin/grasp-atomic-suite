module load_mcp_records_into_memory_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !   Written by  G. Gaigalas               Vilnius, September 2021
        subroutine load_mcp_records_into_memory(nfile, total_secondary_contribution_count, lab, ncontr)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: int64
            implicit none (type, external)

            ! Arguments
            integer(int64), intent(in) :: total_secondary_contribution_count
            integer, intent(in)       :: nfile
            integer, intent(out)      :: lab
            integer, intent(out)      :: ncontr
        end subroutine load_mcp_records_into_memory
    end interface
end module load_mcp_records_into_memory_interface
