module clear_radial_inhomogeneous_terms_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine clear_radial_inhomogeneous_terms (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer :: orbital_index
        end subroutine clear_radial_inhomogeneous_terms
    end interface
end module clear_radial_inhomogeneous_terms_interface
