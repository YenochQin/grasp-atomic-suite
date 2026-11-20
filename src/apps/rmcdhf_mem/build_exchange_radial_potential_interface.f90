module build_exchange_radial_potential_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_exchange_radial_potential (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            !...This routine performs I/O.
        end subroutine build_exchange_radial_potential
    end interface
end module build_exchange_radial_potential_interface
