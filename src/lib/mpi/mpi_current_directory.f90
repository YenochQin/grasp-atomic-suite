module mpi_current_directory
    implicit none (type, external)

    ! Module interface
    private
    public :: sys_getwd
contains
    ! Return the working directory and its trimmed length. Report errors through stopmpi.
    subroutine sys_getwd (directory_path,path_length)

        ! Module dependencies
        use mpi_state, only: mpi_rank
        use posix_services, only: working_directory

        use mpi_runtime, only: stopmpi
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(out):: directory_path
        integer, intent(out):: path_length

        ! Local declarations
        integer :: getcwd_status

        getcwd_status = working_directory(directory_path)
        path_length = len_trim(directory_path)

        if (getcwd_status/=0) then
            call stopmpi ("sys_getwd: status.ne.0; myid = ", mpi_rank)

        end if
        return
    end subroutine sys_getwd
end module mpi_current_directory
