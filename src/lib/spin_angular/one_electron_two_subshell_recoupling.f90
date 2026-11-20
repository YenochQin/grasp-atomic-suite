! Recouple a one-electron tensor connecting two active subshells.
module one_electron_two_subshell_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: recouple_one_electron_two_subshells
contains
    subroutine recouple_one_electron_two_subshells(active_shell_count,first_shell_position,second_shell_position,two_rank_1, &
        & two_rank_2, &
        tensor_rank,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: one
        use csf_pair_angular_state, only: peel_orbital_indices, bra_subshell_quantum_codes, peel_orbital_count
        use first_active_shell_recoupling, only: recouple_first_active_shell
        use spectator_shell_chain_recoupling, only: recouple_spectator_shell_chain
        use pair_tensor_recoupling, only: recouple_pair_tensor
        use final_shell_recoupling, only: recouple_final_shell
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: active_shell_count, first_shell_position, second_shell_position, two_rank_1, &
            two_rank_2, tensor_rank, evaluate_coefficient
        integer, intent(out)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer :: orbital_index, coupling_step_count, coupled_shell_count
        real(real64) :: first_bra_j_weight, bra_j_weight_product, diagram_factor
        coupling_allowed=1
        orbital_index=peel_orbital_indices(first_shell_position)
        first_bra_j_weight=dble(bra_subshell_quantum_codes(3,orbital_index))
        orbital_index=peel_orbital_indices(second_shell_position)
        bra_j_weight_product=dble(bra_subshell_quantum_codes(3,orbital_index))
        bra_j_weight_product=first_bra_j_weight*bra_j_weight_product
        recoupling_coefficient=one/dsqrt(bra_j_weight_product)
        if(peel_orbital_count == 1 .and. active_shell_count == -1)return
        ! A count of -1 extends recoupling through the entire peel chain.
        if(active_shell_count == -1) then
            coupled_shell_count = peel_orbital_count
        else
            coupled_shell_count = active_shell_count
        end if
        coupling_allowed=0
        coupling_step_count=coupled_shell_count-second_shell_position
        if(coupling_step_count > 1) then
            call recouple_spectator_shell_chain(second_shell_position,coupled_shell_count,2*tensor_rank,evaluate_coefficient, &
                & coupling_allowed, &
                diagram_factor)
            if(coupling_allowed == 0)return
            recoupling_coefficient=diagram_factor*recoupling_coefficient
            coupling_allowed=0
        end if
        if(second_shell_position /= coupled_shell_count) then
            call recouple_final_shell(coupled_shell_count,second_shell_position,2*tensor_rank,evaluate_coefficient, &
                & coupling_allowed, &
                diagram_factor)
            if(coupling_allowed == 0)return
            recoupling_coefficient=diagram_factor*recoupling_coefficient
            coupling_allowed=0
        end if
        ! Beyond the second active shell the tensor carries its combined rank. Between the
        ! active shells it carries only the first operator's rank.
        call recouple_pair_tensor(first_shell_position,second_shell_position,two_rank_1,two_rank_2,2*tensor_rank, &
            evaluate_coefficient,coupling_allowed,diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        if(first_shell_position == 1.and.second_shell_position == 2)return
        coupling_allowed=0
        call recouple_first_active_shell(first_shell_position,two_rank_1,evaluate_coefficient,coupling_allowed,diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        coupling_step_count=second_shell_position-first_shell_position
        if(first_shell_position == 1)coupling_step_count=second_shell_position-1-first_shell_position
        if(coupling_step_count <= 1)return
        coupling_allowed=0
        call recouple_spectator_shell_chain(first_shell_position,second_shell_position,two_rank_1,evaluate_coefficient, &
            & coupling_allowed, &
            diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        return
    end subroutine recouple_one_electron_two_subshells
end module one_electron_two_subshell_recoupling
