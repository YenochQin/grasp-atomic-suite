! Open and rewind a formatted or unformatted file and return an open status.
module file_open
    implicit none (type, external)

    ! Module interface
    private
    public :: open_calculation_file
contains
    subroutine open_calculation_file(file_unit, file_path, record_format, requested_status, open_status)
        ! Module dependencies
        use standard_io_units, only: stderr_unit
        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: file_unit
        integer, intent(out) :: open_status
        character (len = *), intent(in) :: file_path
        character (len = *), intent(in) :: record_format
        character (len = *), intent(in) :: requested_status

        ! Local declarations
        integer :: io_status, path_length
        ! Open with status UNKNOWN. requested_status is retained for error diagnostics.
        file_unit = -1
        open(newunit=file_unit, file=file_path, form=record_format, status="UNKNOWN", iostat=io_status, &
            position="asis")
        if (io_status == 0) then
            rewind (file_unit)
            open_status = 0
        else
            path_length = len_trim(file_path)
            write (stderr_unit, *) "OPENFL: Error opening file ", file_path(1:path_length), " as "&
                , requested_status, ";"
            write (stderr_unit, *) "The argument RSTAT=", requested_status, " is not used !"
            open_status = 1
        end if
        return
    end subroutine open_calculation_file
end module file_open
