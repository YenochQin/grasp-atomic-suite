! Test whether a character is one of the ten decimal digits.
module decimal_digit_check
    implicit none (type, external)

    ! Module interface
    private
    public :: is_decimal_digit
contains
    logical function is_decimal_digit (character_to_test)
        implicit none (type, external)

        ! Arguments
        character , intent(in) :: character_to_test

        ! Local declarations
        integer :: digit_value
        character :: decimal_digits(0:9)

        ! Initial values
        data decimal_digits/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"/
        do digit_value = 0, 9
            if (character_to_test /= decimal_digits(digit_value)) cycle
            is_decimal_digit = .true.
            return
        end do
        is_decimal_digit = .false.
    end function is_decimal_digit
end module decimal_digit_check
