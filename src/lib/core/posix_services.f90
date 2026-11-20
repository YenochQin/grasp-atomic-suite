! Wrap POSIX process, filesystem, and directory operations for Fortran callers.
module posix_services
    use iso_c_binding, only: c_char, c_int, c_size_t, c_ptr, c_null_char, c_associated
    implicit none (type, external)
    private
    public :: run_shell, change_directory, working_directory, path_access, make_directory
    public :: terminate_process

    interface
        subroutine c_exit(status) bind(c, name="exit")
            import c_int
            implicit none (type, external)
            integer(c_int), value :: status
        end subroutine c_exit
        function c_system(command) bind(c, name="system") result(status)
            import c_char, c_int
            implicit none (type, external)
            character(kind=c_char), intent(in) :: command(*)
            integer(c_int) :: status
        end function c_system
        function c_chdir(path) bind(c, name="chdir") result(status)
            import c_char, c_int
            implicit none (type, external)
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
        function c_getcwd(buffer, capacity) bind(c, name="getcwd") result(address)
            import c_char, c_size_t, c_ptr
            implicit none (type, external)
            character(kind=c_char), intent(out) :: buffer(*)
            integer(c_size_t), value :: capacity
            type(c_ptr) :: address
        end function c_getcwd
        function c_access(path, mode) bind(c, name="access") result(status)
            import c_char, c_int
            implicit none (type, external)
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value :: mode
            integer(c_int) :: status
        end function c_access
    end interface
contains
    subroutine terminate_process(status)
        integer, intent(in), optional :: status
        if (present(status)) then
            call c_exit(int(status,c_int))
        else
            call c_exit(0_c_int)
        end if
    end subroutine terminate_process

    ! Return the raw wait status from POSIX system.
    integer function run_shell(command) result(status)
        character(*), intent(in) :: command
        character(kind=c_char) :: c_command(len_trim(command)+1)
        integer :: character_index, command_length
        command_length = len_trim(command)
        do character_index = 1, command_length
            c_command(character_index) = command(character_index:character_index)
        end do
        c_command(command_length+1) = c_null_char
        status = c_system(c_command)
    end function run_shell

    integer function change_directory(path) result(status)
        character(*), intent(in) :: path
        character(kind=c_char) :: c_path(len_trim(path)+1)
        integer :: character_index, path_length
        path_length = len_trim(path)
        do character_index = 1, path_length
            c_path(character_index) = path(character_index:character_index)
        end do
        c_path(path_length+1) = c_null_char
        status = c_chdir(c_path)
    end function change_directory

    integer function working_directory(path) result(status)
        character(*), intent(out) :: path
        character(kind=c_char) :: c_path(len(path)+1)
        type(c_ptr) :: buffer_address
        integer :: character_index
        path = ""
        buffer_address = c_getcwd(c_path, int(size(c_path),c_size_t))
        status = 1
        if (.not.c_associated(buffer_address)) return
        do character_index = 1, len(path)
            if (c_path(character_index) == c_null_char) exit
            path(character_index:character_index) = c_path(character_index)
        end do
        status = 0
    end function working_directory

    integer function path_access(path, mode) result(status)
        character(*), intent(in) :: path, mode
        character(kind=c_char) :: c_path(len_trim(path)+1)
        integer(c_int) :: access_mask
        integer :: character_index, path_length
        access_mask = 0_c_int
        do character_index = 1, len_trim(mode)
            select case(mode(character_index:character_index))
            case("r")
                access_mask = ior(access_mask,4_c_int)
            case("w")
                access_mask = ior(access_mask,2_c_int)
            case("x")
                access_mask = ior(access_mask,1_c_int)
            case(" ")
            case default
                status = -1
                return
            end select
        end do
        path_length = len_trim(path)
        do character_index = 1, path_length
            c_path(character_index) = path(character_index:character_index)
        end do
        c_path(path_length+1) = c_null_char
        status = c_access(c_path,access_mask)
    end function path_access

    ! Quote the path as one shell argument, including any literal single quotes.
    subroutine make_directory(path, status)
        character(*), intent(in) :: path
        integer, intent(out) :: status
        character(*), parameter :: command_prefix = "mkdir -p -m 775 -- '"
        character(5*len_trim(path)+len(command_prefix)+1) :: command
        integer :: character_index, command_length
        command = command_prefix
        command_length = len(command_prefix)
        do character_index = 1, len_trim(path)
            if (path(character_index:character_index) == "'") then
                command(command_length+1:command_length+5) = "'"//achar(34)//"'"//achar(34)//"'"
                command_length = command_length+5
            else
                command(command_length+1:command_length+1) = path(character_index:character_index)
                command_length = command_length+1
            end if
        end do
        command(command_length+1:command_length+1) = "'"
        status = run_shell(command(:command_length+1))
    end subroutine make_directory
end module posix_services
