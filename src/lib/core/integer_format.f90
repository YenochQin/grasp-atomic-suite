! Convert an integer to text and return the width of the formatted value.
module integer_format
    implicit none (type, external)
    private
    public :: convrt
contains
      subroutine convrt(integer_value, number_text, text_length)
      implicit none (type, external)
      integer,   intent(in)    :: integer_value
      integer,   intent(out)   :: text_length
      character, intent(inout) :: number_text*(*)
      character(len=6)                  :: integer_format
      character(len=2) :: two_digit_widths(0:10)
      character           :: one_digit_widths(9)
      data one_digit_widths /"1","2","3","4","5","6","7","8","9"/
      data two_digit_widths /"10","11","12","13","14","15","16","17","18","19","20"/
      if (integer_value < 0) then
         text_length = log10(dble((-integer_value))) + 2
      else if (integer_value == 0) then
         text_length = 1
      else
         text_length = log10(dble(integer_value)) + 1
      end if
      if (text_length > len(number_text)) then
         write (6, *) "CONVRT: Length of CNUM inadeuate."
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
      end subroutine convrt
end module integer_format
