subroutine merge_transition_basis_records(state_name,initial_core_subshell_count,final_core_subshell_count)
    ! Module dependencies
    use csf_basis_merging, only: merge_csf_bases
    implicit none (type, external)

    ! Arguments
    character(len=24), intent(in) :: state_name(2)
    integer, intent(in) :: initial_core_subshell_count,final_core_subshell_count

    ! Local declarations
    integer :: status
    character(len=256) :: message
    call merge_csf_bases(initial_core_subshell_count,final_core_subshell_count,status,message)
    if (status /= 0) write(*,*) "CSF descriptor merge ",trim(state_name(1))," / ",trim(state_name(2)),": ",trim(message)
    if (status /= 0) error stop "CSF descriptor merge failed"
end subroutine merge_transition_basis_records
