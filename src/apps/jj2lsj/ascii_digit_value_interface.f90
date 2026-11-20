module ascii_digit_value_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function ascii_digit_value (cst)
            implicit none (type, external)

            ! Arguments
            character (len = 1), intent(in) :: cst
        end function ascii_digit_value
    end interface
end module ascii_digit_value_interface
