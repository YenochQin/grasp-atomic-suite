module orbital_letter_to_l_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function orbital_letter_to_l (symbol)
            implicit none (type, external)

            ! Arguments
            character (len = 1), intent(in) :: symbol
        end function orbital_letter_to_l
    end interface
end module orbital_letter_to_l_interface
