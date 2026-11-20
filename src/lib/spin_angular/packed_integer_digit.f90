! Extract a digit from a nonnegative packed integer in the specified radix.
module packed_integer_digit
    implicit none (type, external)
    private
    public :: jthn
contains
      integer function jthn(packed_value,digit_position,radix) result(digit_value)
      implicit none (type, external)
      integer, intent(in) :: radix, digit_position, packed_value
      if(digit_position == 1) then
        digit_value=mod(packed_value,radix)
      else if(digit_position == 2) then
        digit_value=mod(packed_value/radix,radix)
      else if(digit_position == 3) then
        digit_value=mod(packed_value/(radix*radix),radix)
      else if(digit_position == 4) then
        digit_value=mod(packed_value/(radix*radix*radix),radix)
      else if(digit_position == 5) then
        digit_value=mod(packed_value/(radix*radix*radix*radix),radix)
      else
        write(6,"(A)") " ERROR IN JTHN "
        stop
      end if
      return
      end function jthn
end module packed_integer_digit
