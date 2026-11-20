! Capture the system clock count used to time a run.
module run_timer_start
    implicit none (type, external)

    ! Module interface
    private
    public :: begin_run_timer
contains
    subroutine begin_run_timer(start_clock_count, program_name)
        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: start_clock_count
        character (len = *), intent(in) :: program_name

        ! Local declarations
        integer :: clock_tick_rate, maximum_clock_count

        call system_clock (start_clock_count, clock_tick_rate, maximum_clock_count)

        return
    end subroutine begin_run_timer
end module run_timer_start
