module read_transition_mixing_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_transition_mixing_coefficients (name, inpci, init)
            implicit none (type, external)

            ! Arguments
            character (len = 24) :: name(2)
            integer :: inpci
            integer, intent(in) :: init
        end subroutine read_transition_mixing_coefficients
    end interface
end module read_transition_mixing_coefficients_interface
