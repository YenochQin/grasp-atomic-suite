!***********************************************************************
subroutine open_summary_file(file_path)
    ! Module dependencies
    use calculation_file_units, only: summary_output_unit
    !***********************************************************************
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    character  :: file_path*(*)

    ! Local constants
    character(len=9), parameter :: form = "FORMATTED"
    character(len=3), parameter :: status = "NEW"

    ! Local declarations
    integer :: open_status
    !-----------------------------------------------

    call openfl (summary_output_unit, file_path, form, status, open_status)
    if (open_status /= 0) then
        write (stderr_unit, *) "Error when opening ", file_path
        error stop
    end if

    return
end subroutine open_summary_file
