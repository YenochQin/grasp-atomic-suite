module read_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_initial_radial_wavefunctions (state_name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: state_name
        end subroutine read_initial_radial_wavefunctions
    end interface
end module read_initial_radial_wavefunctions_interface
