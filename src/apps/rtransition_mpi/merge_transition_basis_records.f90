module merge_transition_basis_records_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: merge_transition_basis_records
contains

    subroutine merge_transition_basis_records(state_name,initial_core_subshell_count,final_core_subshell_count)
        ! Module dependencies
        use csf_basis_merging, only: merge_csf_bases
        use mpi_csf_collectives, only: require_collective_success
        implicit none (type, external)

        ! Arguments
        character(len=128), intent(in) :: state_name(2)
        integer, intent(in) :: initial_core_subshell_count,final_core_subshell_count

        ! Local declarations
        integer :: status
        character(len=256) :: message
        call merge_csf_bases(initial_core_subshell_count,final_core_subshell_count,status,message)
        if (status /= 0) write(*,*) "CSF descriptor merge ",trim(state_name(1))," / ",trim(state_name(2)),": ",trim(message)
        call require_collective_success(status)
    end subroutine merge_transition_basis_records
end module merge_transition_basis_records_interface
