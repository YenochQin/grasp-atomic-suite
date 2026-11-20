module write_transition_angular_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_transition_angular_coefficients (operator_parity, transition_type_index, name, transition_angular_unit)
            implicit none (type, external)

            ! Arguments
            integer, intent(out) :: transition_angular_unit
            integer, intent(out) :: operator_parity
            integer, intent(in) :: transition_type_index
            character (len = 24) :: name(2)
        end subroutine write_transition_angular_coefficients
    end interface
end module write_transition_angular_coefficients_interface
