module write_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_initial_radial_wavefunctions(file_path)
            implicit none (type, external)

            ! Arguments
            character(len=*), intent(in) :: file_path
        end subroutine write_initial_radial_wavefunctions
    end interface
end module write_initial_radial_wavefunctions_interface
