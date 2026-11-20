module build_direct_radial_potential_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_direct_radial_potential (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
        end subroutine build_direct_radial_potential
    end interface
end module build_direct_radial_potential_interface
