module open_summary_file_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_summary_file
contains

    !***********************************************************************
    subroutine open_summary_file(file_path)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        !***********************************************************************
        use standard_io_units, only: stderr_unit
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        character, intent(in) :: file_path*(*)

        ! Local constants
        character(len=9), parameter :: summary_file_format = "FORMATTED"
        character(len=3), parameter :: summary_file_status = "NEW"

        ! Local declarations
        integer :: open_status
        !-----------------------------------------------

        call open_calculation_file (summary_output_unit, file_path, summary_file_format, summary_file_status, open_status)
        if (open_status /= 0) then
            write (stderr_unit, *) "Error when opening ", file_path
            error stop
        end if

        return
    end subroutine open_summary_file
end module open_summary_file_interface
