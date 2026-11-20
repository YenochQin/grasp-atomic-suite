! Return the current local date and time as fixed-width strings from DATE_AND_TIME.
module calendar_clock
    implicit none (type, external)

    ! Module interface
    private
    public :: read_calendar_clock
contains
    subroutine read_calendar_clock(time_text, date_text)
        implicit none (type, external)

        ! Arguments
        character(len=10), intent(out) :: time_text
        character(len=8), intent(out) :: date_text

        call date_and_time(date_text, time_text)

        return
    end subroutine read_calendar_clock
end module calendar_clock
