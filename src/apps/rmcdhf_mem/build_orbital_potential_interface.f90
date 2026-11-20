module build_orbital_potential_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_orbital_potential (extended_optimal_level, orbital_index, radial_point_count)
            implicit none (type, external)

            ! Arguments
            logical :: extended_optimal_level
            integer :: orbital_index
            integer :: radial_point_count
        end subroutine build_orbital_potential
    end interface
end module build_orbital_potential_interface
