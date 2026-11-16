module input_io
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none (type, external)
    private
    public :: fail_input_io

    abstract interface
        subroutine abort_input()
            implicit none (type, external)
        end subroutine abort_input
    end interface

contains

    subroutine fail_input_io(status, unit, context, abort_handler)
        integer, intent(in) :: status, unit
        character(*), intent(in) :: context
        procedure(abort_input), optional :: abort_handler
        integer :: inquiry_status
        logical :: opened, named
        character(1024) :: filename

        ! Zero denotes an existing input consistency check, not a failed READ.
        if (status == 0) then
            write(error_unit, "(a,a,a,i0,a)") "Input I/O failure: ", trim(context), &
                "; unit=", unit, "; input validation failure"
        else
            write(error_unit, "(a,a,a,i0,a,i0)") "Input I/O failure: ", trim(context), &
                "; unit=", unit, "; iostat=", status
        end if
        inquire(unit=unit, opened=opened, named=named, name=filename, iostat=inquiry_status)
        if (inquiry_status == 0) then
            if (opened) then
                if (named) write(error_unit, "(a,a)") "Input file: ", trim(filename)
            end if
        end if
        flush(error_unit)
        if (present(abort_handler)) call abort_handler()
        error stop 1
    end subroutine fail_input_io

end module input_io
