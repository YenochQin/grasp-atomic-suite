module estimate_from_saved_orbitals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine estimate_from_saved_orbitals (matched_orbital_indices, matched_subshell_count, file_path)
            ! Module dependencies
            use grasp_parameters, only: radial_workspace_capacity, max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: matched_orbital_indices(max_orbitals)
            integer, intent(in) :: matched_subshell_count
            character (len = *), intent(in) :: file_path
        end subroutine estimate_from_saved_orbitals
    end interface
end module estimate_from_saved_orbitals_interface
