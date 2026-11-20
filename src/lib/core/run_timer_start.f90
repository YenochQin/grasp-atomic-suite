! Capture the system clock count used to time a run.
module run_timer_start
    implicit none (type, external)

    ! Module interface
    private
    public :: starttime
contains
    subroutine starttime(start_clock_count, program_name)
        implicit none (type, external)

        ! Arguments
        integer, intent(out) :: start_clock_count
        character (len = *), intent(in) :: program_name

        ! Local declarations
        integer :: nymduhmsm(8)
        integer :: clock_tick_rate, maximum_clock_count
        character :: chdate*8, chtime*10, chzone*5, msg*80

        call system_clock (start_clock_count, clock_tick_rate, maximum_clock_count)

        return
    end subroutine starttime
end module run_timer_start
