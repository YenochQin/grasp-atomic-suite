module mpi_make_directory
    implicit none (type, external)

    ! Module interface
    private
    public :: create_mpi_directory
contains
    ! Create a directory from the first path_length characters and return the system status.
    subroutine create_mpi_directory (directory_path, path_length, directory_status)
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
            print*, " create_mpi_directory error message = ", trim(mkdir_error_message)
        end if
        return
    end subroutine create_mpi_directory
end module mpi_make_directory
