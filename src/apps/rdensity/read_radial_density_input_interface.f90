module read_radial_density_input_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_radial_density_input (name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name
        end subroutine read_radial_density_input
    end interface
end module read_radial_density_input_interface
