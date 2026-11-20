module load_final_transition_basis_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine load_final_transition_basis (final_core_subshell_count, state_name)
            implicit none (type, external)

            ! Arguments
            integer, intent(out) :: final_core_subshell_count
            character (len = 128), intent(in) :: state_name
        end subroutine load_final_transition_basis
    end interface
end module load_final_transition_basis_interface
