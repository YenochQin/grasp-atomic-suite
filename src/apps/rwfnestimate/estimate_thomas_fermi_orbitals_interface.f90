module estimate_thomas_fermi_orbitals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine estimate_thomas_fermi_orbitals (matched_orbital_indices, matched_subshell_count)
            ! Module dependencies
            use grasp_parameters, only: max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: matched_orbital_indices(max_orbitals)
            integer, intent(in) :: matched_subshell_count
        end subroutine estimate_thomas_fermi_orbitals
    end interface
end module estimate_thomas_fermi_orbitals_interface
