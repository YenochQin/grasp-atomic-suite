module read_state_mixing_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_state_mixing_coefficients (name, inpci, iblk)
            implicit none (type, external)

            ! Arguments
            character(len=*), intent(in) :: name
            integer, intent(in) :: inpci
            integer, intent(in) :: iblk
        end subroutine read_state_mixing_coefficients
    end interface
end module read_state_mixing_coefficients_interface
