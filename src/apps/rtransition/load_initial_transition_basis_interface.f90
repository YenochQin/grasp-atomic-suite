module load_initial_transition_basis_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine load_initial_transition_basis (initial_core_subshell_count,state_name)
            implicit none (type, external)

            ! Arguments
            integer :: initial_core_subshell_count
            character(len=24) :: state_name
        end subroutine load_initial_transition_basis
    end interface
end module load_initial_transition_basis_interface
