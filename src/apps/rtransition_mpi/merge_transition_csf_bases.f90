! Align the initial and final descriptor bases in memory for transition matrix elements.
module merge_transition_csf_bases_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: merge_transition_csf_bases
contains

    subroutine merge_transition_csf_bases(state_name)
        ! Module dependencies

        ! Procedure interfaces
        use load_initial_transition_basis_interface, only: load_initial_transition_basis
        use load_final_transition_basis_interface, only: load_final_transition_basis
        use merge_transition_basis_records_interface, only: merge_transition_basis_records
        implicit none (type, external)

        ! Arguments
        character  :: state_name(2)*128

        ! Local declarations
        integer :: initial_core_subshell_count, final_core_subshell_count

        ! Load the initial-state CSF descriptors.
        call load_initial_transition_basis (initial_core_subshell_count, state_name(1))
        ! Load the final-state CSF descriptors.
        call load_final_transition_basis (final_core_subshell_count, state_name(2))
        ! Merge the two orbital orders and remap the CSF records in memory.
        call merge_transition_basis_records (state_name, initial_core_subshell_count, final_core_subshell_count)

        return
    end subroutine merge_transition_csf_bases
end module merge_transition_csf_bases_interface
