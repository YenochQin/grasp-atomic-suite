! Screen the total J, inactive subshells and intermediate couplings for a one-electron tensor.
module one_electron_tensor_recoupling_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: select_one_electron_tensor_recoupling
contains
    subroutine select_one_electron_tensor_recoupling(active_shell_count,first_shell_position,second_shell_position,tensor_rank, &
        & coupling_allowed)
        ! Module dependencies
        use csf_pair_angular_state,             only: peel_orbital_count, bra_intermediate_couplings, &
            ket_intermediate_couplings, peel_orbital_indices, bra_subshell_quantum_codes, ket_subshell_quantum_codes
        use angular_triangle_and_parity, only: satisfies_angular_triangle_and_parity
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: active_shell_count, first_shell_position, second_shell_position, tensor_rank
        integer, intent(out) :: coupling_allowed

        ! Local declarations
        integer :: shell_position, coupling_step_count, orbital_index, prefix_shell_position, suffix_step, &
            prefix_coupling_position, coupled_shell_count, suffix_coupling_position, two_bra_chain_j, two_ket_chain_j
        coupling_allowed=1
        if(peel_orbital_count == 1 .and. active_shell_count == -1)return
        coupling_allowed=0
        if(active_shell_count == -1) then
            coupled_shell_count = peel_orbital_count
        else
            coupled_shell_count = active_shell_count
        end if
        if(coupled_shell_count > 1) then
            two_bra_chain_j=bra_intermediate_couplings(coupled_shell_count-1)-1
            two_ket_chain_j=ket_intermediate_couplings(coupled_shell_count-1)-1
        else if(coupled_shell_count == 1) then
            two_bra_chain_j=bra_intermediate_couplings(coupled_shell_count)-1
            two_ket_chain_j=ket_intermediate_couplings(coupled_shell_count)-1
        else
            print*, "ERROR in RECOP00"
            stop
        end if

        ! Reject total-J couplings that cannot accommodate the tensor rank.
        if(satisfies_angular_triangle_and_parity(two_bra_chain_j,two_ket_chain_j,2*tensor_rank) == 0)return
        coupling_allowed=1
        if(peel_orbital_count == 1)return
        do shell_position=1,peel_orbital_count
            if(first_shell_position /= shell_position) then
                if(second_shell_position /= shell_position) then
                    orbital_index=peel_orbital_indices(shell_position)
                    if(bra_subshell_quantum_codes(1,orbital_index) /= ket_subshell_quantum_codes(1, &
                        orbital_index))coupling_allowed=0
                    if(bra_subshell_quantum_codes(2,orbital_index) /= ket_subshell_quantum_codes(2, &
                        orbital_index))coupling_allowed=0
                    if(bra_subshell_quantum_codes(3,orbital_index) /= ket_subshell_quantum_codes(3, &
                        orbital_index))coupling_allowed=0
                end if
            end if
        end do
        if(coupling_allowed == 0)return
        if(coupled_shell_count <= 2)return
        if(first_shell_position <= 2)return
        do prefix_shell_position=3,first_shell_position
            prefix_coupling_position=prefix_shell_position-2
            if(bra_intermediate_couplings(prefix_coupling_position) /= ket_intermediate_couplings( &
                prefix_coupling_position))coupling_allowed=0
            if(coupling_allowed == 0)return
        end do
        coupling_step_count=coupled_shell_count-second_shell_position
        if(coupling_step_count > 0) then
            do suffix_step=1,coupling_step_count
                suffix_coupling_position=second_shell_position-2+suffix_step
                two_bra_chain_j=bra_intermediate_couplings(suffix_coupling_position)-1
                two_ket_chain_j=ket_intermediate_couplings(suffix_coupling_position)-1
                if(satisfies_angular_triangle_and_parity(two_bra_chain_j,two_ket_chain_j,2*tensor_rank) == 0)coupling_allowed=0
                if(coupling_allowed == 0)return
            end do
        end if
        return
    end subroutine select_one_electron_tensor_recoupling
end module one_electron_tensor_recoupling_selection
