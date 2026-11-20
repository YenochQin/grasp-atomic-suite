module estimate_orbital_energy_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine estimate_orbital_energy (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
        end subroutine estimate_orbital_energy
    end interface
end module estimate_orbital_energy_interface
