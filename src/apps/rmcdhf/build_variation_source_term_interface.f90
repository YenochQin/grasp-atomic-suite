module build_variation_source_term_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_variation_source_term (orbital_index)
            implicit none (type, external)

            ! Arguments
            integer :: orbital_index
        end subroutine build_variation_source_term
    end interface
end module build_variation_source_term_interface
