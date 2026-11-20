! Convert twice an angular momentum to integer or numerator/2 text.
module half_integer_format
    implicit none (type, external)

    ! Module interface
    private
    public :: format_half_integer
contains
    subroutine format_half_integer(twice_angular_momentum, angular_momentum_text, text_length)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit

        implicit none (type, external)

        ! Arguments
        integer,   intent(in)    :: twice_angular_momentum
        integer,   intent(out)   :: text_length
        character, intent(inout) :: angular_momentum_text*(*)

        ! Local declarations
        integer                           :: numerator_to_print
        character(len=6)                  :: integer_format
        character(len=2) :: two_digit_widths(0:10)
        character           :: one_digit_widths(9)

        ! Initial values
        data one_digit_widths /"1","2","3","4","5","6","7","8","9"/
        data two_digit_widths /"10","11","12","13","14","15","16","17","18","19","20"/
        if(mod(twice_angular_momentum,2) == 0) then
            numerator_to_print = twice_angular_momentum/2
        else
            numerator_to_print = twice_angular_momentum
        end if
        if (numerator_to_print < 0) then
            text_length = log10(dble((-numerator_to_print))) + 2
        else if (numerator_to_print == 0) then
            text_length = 1
        else
            text_length = log10(dble(numerator_to_print)) + 1
        end if
        if (text_length > len(angular_momentum_text)) then
            write (output_unit, *) "CONVRT_DOUBLE: Length of CNUM inadeuate."
            error stop
        else
            if (text_length <= 9) then
                integer_format = "(1I"//one_digit_widths(text_length)//")"
                write (angular_momentum_text(1:text_length), integer_format(1:5)) numerator_to_print
            else
                integer_format = "(1I"//two_digit_widths(text_length-10)//")"
                write (angular_momentum_text(1:text_length), integer_format(1:6)) numerator_to_print
            end if
            if(mod(twice_angular_momentum,2) /= 0) then
                if (text_length+2 > len(angular_momentum_text)) then
                    write (output_unit, *) "CONVRT_DOUBLE: Length of CNUM inadeuate."
                    error stop
                else
                    angular_momentum_text(1:text_length+2) = angular_momentum_text(1:text_length)//"/2"
                    text_length = text_length + 2
                end if
            end if
        end if
        return
    end subroutine format_half_integer
end module half_integer_format
