module mpi_make_directory
    implicit none (type, external)

    ! Module interface
    private
    public :: sys_mkdir
contains
    ! Create a directory from the first path_length characters and return the system status.
    subroutine sys_mkdir (directory_path, path_length, directory_status)
        ! Module dependencies
        use posix_services, only: make_directory
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in):: directory_path
        integer, intent(in):: path_length
        integer, intent(out):: directory_status

        ! Local declarations
        character(len=255) :: mkdir_error_message

        mkdir_error_message = " no message passed from the shell "
        directory_status = -999

        call make_directory (directory_path(1:path_length), directory_status)
        if (directory_status/=0) then
            print*, " sys_mkdir error message = ", trim(mkdir_error_message)
        end if
        return
    end subroutine sys_mkdir
end module mpi_make_directory
