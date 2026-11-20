module merge_transition_basis_records_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine merge_transition_basis_records (state_name, initial_core_subshell_count, final_core_subshell_count)
            implicit none (type, external)

            ! Arguments
            character (len = 24), intent(in) :: state_name(2)
            integer, intent(in) :: initial_core_subshell_count
            integer, intent(in) :: final_core_subshell_count
        end subroutine merge_transition_basis_records
    end interface
end module merge_transition_basis_records_interface
