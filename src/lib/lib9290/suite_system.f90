! Standard Fortran interfaces to the POSIX services used by the suite.
! Keep SYSTEM's raw wait status (e.g. exit 7 -> 1792), rather than silently
! changing it to EXECUTE_COMMAND_LINE's unencoded exit status.
module suite_system
    use iso_c_binding, only: c_char, c_int, c_size_t, c_ptr, c_null_char, c_associated
    implicit none
    private
    public :: run_shell, change_directory, working_directory, path_access, make_directory
    public :: terminate_process

    interface
        subroutine c_exit(status) bind(C, name='exit')
            import c_int
            integer(c_int), value :: status
        end subroutine c_exit
        function c_system(command) bind(C, name='system') result(status)
            import c_char, c_int
            character(kind=c_char), intent(in) :: command(*)
            integer(c_int) :: status
        end function c_system
        function c_chdir(path) bind(C, name='chdir') result(status)
            import c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: status
        end function c_chdir
        function c_getcwd(buffer, capacity) bind(C, name='getcwd') result(address)
            import c_char, c_size_t, c_ptr
            character(kind=c_char), intent(out) :: buffer(*)
            integer(c_size_t), value :: capacity
            type(c_ptr) :: address
        end function c_getcwd
        function c_access(path, mode) bind(C, name='access') result(status)
            import c_char, c_int
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

    integer function run_shell(command) result(status)
        character(*), intent(in) :: command
        character(kind=c_char) :: buffer(len_trim(command)+1)
        integer :: i, n
        n = len_trim(command)
        do i = 1, n
            buffer(i) = command(i:i)
        end do
        buffer(n+1) = c_null_char
        status = c_system(buffer)
    end function run_shell

    integer function change_directory(path) result(status)
        character(*), intent(in) :: path
        character(kind=c_char) :: buffer(len_trim(path)+1)
        integer :: i, n
        n = len_trim(path)
        do i = 1, n
            buffer(i) = path(i:i)
        end do
        buffer(n+1) = c_null_char
        status = c_chdir(buffer)
    end function change_directory

    integer function working_directory(path) result(status)
        character(*), intent(out) :: path
        character(kind=c_char) :: buffer(len(path)+1)
        type(c_ptr) :: address
        integer :: i
        path = ''
        address = c_getcwd(buffer, int(size(buffer),c_size_t))
        status = 1
        if (.not.c_associated(address)) return
        do i = 1, len(path)
            if (buffer(i) == c_null_char) exit
            path(i:i) = buffer(i)
        end do
        status = 0
    end function working_directory

    integer function path_access(path, mode) result(status)
        character(*), intent(in) :: path, mode
        character(kind=c_char) :: buffer(len_trim(path)+1)
        integer(c_int) :: mask
        integer :: i, n
        mask = 0_c_int
        do i = 1, len_trim(mode)
            select case(mode(i:i))
            case('r')
                mask = ior(mask,4_c_int)
            case('w')
                mask = ior(mask,2_c_int)
            case('x')
                mask = ior(mask,1_c_int)
            case(' ')
            case default
                status = -1
                return
            end select
        end do
        n = len_trim(path)
        do i = 1, n
            buffer(i) = path(i:i)
        end do
        buffer(n+1) = c_null_char
        status = c_access(buffer,mask)
    end function path_access

    subroutine make_directory(path, status)
        character(*), intent(in) :: path
        integer, intent(out) :: status
        character(*), parameter :: prefix = "mkdir -p -m 775 -- '"
        character(5*len_trim(path)+len(prefix)+1) :: command
        integer :: i, used
        command = prefix
        used = len(prefix)
        do i = 1, len_trim(path)
            if (path(i:i) == "'") then
                command(used+1:used+5) = "'"//achar(34)//"'"//achar(34)//"'"
                used = used+5
            else
                command(used+1:used+1) = path(i:i)
                used = used+1
            end if
        end do
        command(used+1:used+1) = "'"
        status = run_shell(command(:used+1))
    end subroutine make_directory
end module suite_system
