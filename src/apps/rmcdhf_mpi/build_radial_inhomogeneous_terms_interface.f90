module build_radial_inhomogeneous_terms_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_radial_inhomogeneous_terms (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer :: orbital_index
        end subroutine build_radial_inhomogeneous_terms
    end interface
end module build_radial_inhomogeneous_terms_interface
