module ascii_digit_value_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: ascii_digit_value
contains

    !***********************************************************************
    !                                                                      *
    integer function ascii_digit_value (digit_character)
        !                                                                      *
        !                                                                      *
        !                                                                      *
        !***********************************************************************
        !-----------------------------------------------
        implicit none (type, external)

        ! Arguments
        character , intent(in) :: digit_character

        ! Local declarations
        integer :: digit_value
        character :: cdgt(0:9)

        ! Initial values
        !-----------------------------------------------
        !
        data cdgt/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"/
        !
        do digit_value = 0, 9
            if (digit_character /= cdgt(digit_value)) cycle
            ascii_digit_value = digit_value
            exit
        end do
        !
        return
    end function ascii_digit_value
end module ascii_digit_value_interface
