module write_transformed_radial_wavefunctions_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_transformed_radial_wavefunctions (name)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: name(2)
        end subroutine write_transformed_radial_wavefunctions
    end interface
end module write_transformed_radial_wavefunctions_interface
