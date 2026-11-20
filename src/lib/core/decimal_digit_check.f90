! Test whether a character is one of the ten decimal digits.
module decimal_digit_check
    implicit none (type, external)
    private
    public :: ldigit
contains
      logical function ldigit (character_to_test)
      implicit none (type, external)
      character , intent(in) :: character_to_test
      integer :: digit_value
      character :: decimal_digits(0:9)
      data decimal_digits/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"/
      do digit_value = 0, 9
         if (character_to_test /= decimal_digits(digit_value)) cycle
         ldigit = .true.
         goto 2
      end do
      ldigit = .false.
    2 continue
      return
      end function ldigit
end module decimal_digit_check
