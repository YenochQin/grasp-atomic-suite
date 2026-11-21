! Extract a digit from a nonnegative packed integer in the specified radix.
module packed_integer_digit
    implicit none (type, external)

    ! Module interface
    private
    public :: extract_packed_integer_digit
contains
    integer function extract_packed_integer_digit(packed_value,digit_position,radix) result(digit_value)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        implicit none (type, external)

        ! Arguments
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
            write(output_unit,"(A)") " ERROR IN extract_packed_integer_digit "
            stop
        end if
        return
    end function extract_packed_integer_digit
end module packed_integer_digit
