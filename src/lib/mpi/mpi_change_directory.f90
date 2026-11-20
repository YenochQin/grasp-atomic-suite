module mpi_change_directory
    implicit none (type, external)

    ! Module interface
    private
    public :: change_mpi_directory
contains
    ! Change the working directory using the first path_length characters.
    subroutine change_mpi_directory (directory_path, path_length, directory_status)
        ! Module dependencies
        use posix_services, only: change_directory
        implicit none (type, external)

        ! Arguments
        character(len=*), intent(in):: directory_path
        integer, intent(in):: path_length
        integer, intent(out):: directory_status

        directory_status = change_directory (directory_path(1:path_length))

        return
    end subroutine change_mpi_directory
end module mpi_change_directory
