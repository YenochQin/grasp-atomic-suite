module mpi_current_directory
    implicit none (type, external)

    ! Module interface
    private
    public :: get_mpi_working_directory
contains
    ! Return the working directory and its trimmed length. Report errors through stopmpi.
    subroutine get_mpi_working_directory (directory_path,path_length)

        ! Module dependencies
        use mpi_state, only: mpi_rank
        use posix_services, only: working_directory

        use mpi_runtime, only: finalize_mpi_runtime
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(out):: directory_path
        integer, intent(out):: path_length

        ! Local declarations
        integer :: getcwd_status

        getcwd_status = working_directory(directory_path)
        path_length = len_trim(directory_path)

        if (getcwd_status/=0) then
            call finalize_mpi_runtime ("sys_getwd: status.ne.0; myid = ", mpi_rank)

        end if
        return
    end subroutine get_mpi_working_directory
end module mpi_current_directory
