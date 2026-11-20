! Screen the subshell and intermediate J values for a scalar one-electron term.
module one_electron_scalar_recoupling_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: select_one_electron_scalar_recoupling
contains
    subroutine select_one_electron_scalar_recoupling(active_shell_count,first_shell_position,second_shell_position, &
        & third_shell_position, &
        fourth_shell_position,interaction_case,coupling_allowed)
        ! Module dependencies
        use csf_pair_angular_state,          only: peel_orbital_count, peel_orbital_indices, bra_intermediate_couplings, &
            ket_intermediate_couplings, bra_subshell_quantum_codes, ket_subshell_quantum_codes
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: active_shell_count, first_shell_position, second_shell_position, &
            third_shell_position, fourth_shell_position, interaction_case
        integer, intent(out) :: coupling_allowed

        ! Local declarations
        integer :: shell_position, first_skipped_coupling, coupling_skip_end, orbital_index, orbital_index_1, &
            orbital_index_2, shell_code_slot, coupled_shell_count
        coupling_allowed=1
        if(peel_orbital_count == 1 .and. active_shell_count == -1)return
        if(active_shell_count == -1) then
            coupled_shell_count = peel_orbital_count
        else
            coupled_shell_count = active_shell_count
        end if
        orbital_index_1=peel_orbital_indices(first_shell_position)
        orbital_index_2=peel_orbital_indices(second_shell_position)
        if ((first_shell_position == 1 .and. second_shell_position == 2) .or. interaction_case == 0) then
            do shell_position=1,coupled_shell_count
                orbital_index=peel_orbital_indices(shell_position)
                if(shell_position < coupled_shell_count-1) then
                    if(bra_intermediate_couplings(shell_position) /= ket_intermediate_couplings(shell_position)) &
                        coupling_allowed=0
                end if
                if(interaction_case /= 0) then
                    if(shell_position == first_shell_position) cycle
                    if(shell_position == second_shell_position) cycle
                end if
                do shell_code_slot=1,3
                    if(bra_subshell_quantum_codes(shell_code_slot,orbital_index) /= &
                        ket_subshell_quantum_codes(shell_code_slot,orbital_index))coupling_allowed=0
                end do
            end do
            return
        end if
        do shell_position=1,coupled_shell_count
            orbital_index=peel_orbital_indices(shell_position)
            if(shell_position < coupled_shell_count-1) then
                first_skipped_coupling=first_shell_position-1
                coupling_skip_end=second_shell_position-1
                if(first_shell_position == 1)first_skipped_coupling=first_shell_position
                if (shell_position < first_skipped_coupling .or. shell_position >= coupling_skip_end) then
                    if(bra_intermediate_couplings(shell_position) /= ket_intermediate_couplings(shell_position)) &
                        coupling_allowed=0
                end if
            end if
            if(shell_position == first_shell_position) cycle
            if(shell_position == second_shell_position) cycle
            if((interaction_case == 2).and.(shell_position == third_shell_position)) cycle
            if((interaction_case == 3).and.(shell_position == third_shell_position)) cycle
            if((interaction_case == 3).and.(shell_position == fourth_shell_position)) cycle
            do shell_code_slot=1,3
                if(bra_subshell_quantum_codes(shell_code_slot,orbital_index) /= &
                    ket_subshell_quantum_codes(shell_code_slot,orbital_index))coupling_allowed=0
            end do
        end do
        return
    end subroutine select_one_electron_scalar_recoupling
end module one_electron_scalar_recoupling_selection
