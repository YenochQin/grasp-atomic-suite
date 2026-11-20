module prepare_transition_mixing_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine prepare_transition_mixing_coefficients (name, inpci)
            implicit none (type, external)

            ! Arguments
            character (len = 128), intent(inout) :: name(2)
            integer, intent(in) :: inpci
        end subroutine prepare_transition_mixing_coefficients
    end interface
end module prepare_transition_mixing_coefficients_interface
