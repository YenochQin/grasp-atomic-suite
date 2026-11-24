! Report the elapsed wall-clock time for a serial, OpenMP, or MPI application.
module run_timer
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    implicit none (type, external)

    private
    public :: begin_run_timer, finish_run_timer
contains
    subroutine begin_run_timer (start_clock_count, program_name, log_unit)
        integer(int64), intent(out) :: start_clock_count
        character(len=*), intent(in) :: program_name
        integer, optional, intent(in) :: log_unit

        integer :: unit

        call system_clock (count=start_clock_count)
        unit = output_unit
        if (present(log_unit)) unit = log_unit
        write (unit, *)
        write (unit, "(a)") trim(program_name)
        call write_timestamp (unit, "Date and Time:")
    end subroutine begin_run_timer

    subroutine finish_run_timer (start_clock_count, program_name, log_unit)
        integer(int64), intent(in) :: start_clock_count
        character(len=*), intent(in) :: program_name
        integer, optional, intent(in) :: log_unit

        integer :: unit
        integer(int64) :: finish_clock_count, ticks_per_second, clock_limit, elapsed_ticks

        call system_clock (finish_clock_count, ticks_per_second, clock_limit)
        unit = output_unit
        if (present(log_unit)) unit = log_unit
        write (unit, *)
        call write_timestamp (unit, "Finish Date and Time:")
        if (ticks_per_second > 0 .and. start_clock_count >= 0 .and. finish_clock_count >= 0) then
            if (finish_clock_count >= start_clock_count) then
                elapsed_ticks = finish_clock_count - start_clock_count
            else
                elapsed_ticks = (clock_limit - start_clock_count) + finish_clock_count + 1
            end if
            write (unit, "(a,i0,a)") "Wall time: ", elapsed_ticks / ticks_per_second, " seconds"
        else
            write (unit, "(a)") "Wall time: unavailable"
        end if
        write (unit, "(a)") trim(program_name) // ": Execution complete."
    end subroutine finish_run_timer

    subroutine write_timestamp (unit, label)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: label

        character(len=8) :: date_text
        character(len=10) :: time_text
        character(len=5) :: zone_text

        call date_and_time (date_text, time_text, zone_text)
        write (unit, "(*(a))") label, " ", date_text(:4), "-", date_text(5:6), "-", date_text(7:), &
            " ", time_text(:2), ":", time_text(3:4), ":", time_text(5:), " ", trim(zone_text)
    end subroutine write_timestamp
end module run_timer
