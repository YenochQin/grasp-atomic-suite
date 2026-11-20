! Report elapsed wall-clock seconds and the local completion date and time.
module run_timer_finish
    implicit none (type, external)
    private
    public :: stoptime
contains
      subroutine stoptime(start_clock_count, program_name)
      implicit none (type, external)
      integer , intent(in) :: start_clock_count
      character (len = *), intent(in) :: program_name
      integer :: elapsed_clock_count, clock_tick_rate, maximum_clock_count, elapsed_seconds
      integer :: date_time_values(8)
      character :: date_text*8, time_text*10, time_zone_text*5, elapsed_seconds_text*8, elapsed_time_message*80

      write (6, *)
      write (6, *) "Wall time:"

      call system_clock (elapsed_clock_count, clock_tick_rate, maximum_clock_count)
      elapsed_clock_count = elapsed_clock_count - start_clock_count
      elapsed_seconds = elapsed_clock_count/clock_tick_rate
      write (elapsed_seconds_text, "(I8)") elapsed_seconds
      elapsed_time_message = elapsed_seconds_text//" seconds"
      write (6, *) elapsed_time_message(1:len_trim(elapsed_time_message))

      write (6, *)
      write (6, *) "Finish Date and Time:"

      call date_and_time (date_text, time_text, time_zone_text, date_time_values)

      print*, "  Date (Yr/Mon/Day): ",                     &
              date_text(1:4),"/",date_text(5:6),"/",date_text(7:8)
      print*, "  Time (Hr/Min/Sec): ",                     &
              time_text(1:2),"/",time_text(3:4),"/",time_text(5:10)
      print*, "  Zone: ",time_zone_text

      write (6, *)
      write (6, *) program_name//": Execution complete."

      return
      end subroutine stoptime
end module run_timer_finish
