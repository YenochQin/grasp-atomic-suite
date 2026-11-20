module print_radial_wavefunction_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine print_radial_wavefunction (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
        end subroutine print_radial_wavefunction
    end interface
end module print_radial_wavefunction_interface
