module estimate_hydrogenic_orbitals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine estimate_hydrogenic_orbitals (matched_orbital_indices, matched_subshell_count, modify, nradial)
            ! Module dependencies
            use grasp_parameters, only: max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: matched_orbital_indices(max_orbitals)
            integer, intent(in) :: matched_subshell_count, nradial
            logical, intent(in) :: modify
        end subroutine estimate_hydrogenic_orbitals
    end interface
end module estimate_hydrogenic_orbitals_interface
