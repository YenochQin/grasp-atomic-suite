! Convert an integer to text and include the caller name if the output buffer is too short.
module integer_format_with_context
    implicit none (type, external)
    private
    public :: convrt2
contains
      subroutine convrt2(integer_value, number_text, text_length, caller_name)
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      implicit none (type, external)
      integer, intent(in) :: integer_value
      integer, intent(out) :: text_length
      character, intent(inout) :: number_text*(*)
      character, intent(in) :: caller_name*(*)
      character :: integer_format*6
      character :: two_digit_widths(0:10)*2
      character :: one_digit_widths(9)
      data one_digit_widths/ "1", "2", "3", "4", "5", "6", "7", "8", "9"/
      data two_digit_widths/ "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", &
         "20"/
      if (integer_value < 0) then
         text_length = log10(dble((-integer_value))) + 2
      else if (integer_value == 0) then
         text_length = 1
      else
         text_length = log10(dble(integer_value)) + 1
      end if
      if (text_length > len(number_text)) then
         write (stderr_unit, *) "CONVRT: Length of CNUM inadeuate. (from:", caller_name, ")"
         error stop
      else
         if (text_length <= 9) then
            integer_format = "(1I"//one_digit_widths(text_length)//")"
            write (number_text(1:text_length), integer_format(1:5)) integer_value
         else
            integer_format = "(1I"//two_digit_widths(text_length-10)//")"
            write (number_text(1:text_length), integer_format(1:6)) integer_value
         end if
      end if
      return
      end subroutine convrt2
end module integer_format_with_context
