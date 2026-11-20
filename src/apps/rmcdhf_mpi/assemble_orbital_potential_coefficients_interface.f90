module assemble_orbital_potential_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine assemble_orbital_potential_coefficients (extended_optimal_level, orbital_index)
            implicit none (type, external)

            ! Arguments
            logical :: extended_optimal_level
            integer, intent(in) :: orbital_index
        end subroutine assemble_orbital_potential_coefficients
    end interface
end module assemble_orbital_potential_coefficients_interface
