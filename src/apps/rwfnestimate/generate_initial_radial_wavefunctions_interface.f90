module generate_initial_radial_wavefunctions_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine generate_initial_radial_wavefunctions(default_orbitals)
            implicit none (type, external)

            ! Arguments
            character(len=*), intent(in) :: default_orbitals
        end subroutine generate_initial_radial_wavefunctions
    end interface
end module generate_initial_radial_wavefunctions_interface
