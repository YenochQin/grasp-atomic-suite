! Multiply the coupling-chain factors for three ordered active subshells.
module three_subshell_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: recouple_three_subshells
contains
    subroutine recouple_three_subshells(first_shell_position,second_shell_position,third_shell_position,two_rank_1,two_rank_2, &
        two_tensor_rank,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: one
        use csf_pair_angular_state,             only: peel_orbital_indices, bra_subshell_quantum_codes, ket_subshell_quantum_codes
        use first_active_shell_recoupling, only: recouple_first_active_shell
        use later_active_shell_recoupling, only: recouple_later_active_shell
        use spectator_shell_chain_recoupling, only: recouple_spectator_shell_chain
        use pair_tensor_recoupling, only: recouple_pair_tensor
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: first_shell_position,second_shell_position,third_shell_position,two_rank_1, &
            two_rank_2,two_tensor_rank,evaluate_coefficient
        integer, intent(out)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer      :: two_bra_shell_j_3, two_ket_shell_j_3, orbital_index_1, orbital_index_2, orbital_index_3, &
            coupling_step_count
        real(real64) :: bra_j_weight_product, bra_j_weight_1, bra_j_weight_2, bra_j_weight_3, diagram_factor
        orbital_index_1=peel_orbital_indices(first_shell_position)
        orbital_index_2=peel_orbital_indices(second_shell_position)
        orbital_index_3=peel_orbital_indices(third_shell_position)
        bra_j_weight_1=bra_subshell_quantum_codes(3,orbital_index_1)
        bra_j_weight_2=bra_subshell_quantum_codes(3,orbital_index_2)
        bra_j_weight_3=bra_subshell_quantum_codes(3,orbital_index_3)
        bra_j_weight_product=bra_j_weight_1*bra_j_weight_2*bra_j_weight_3
        recoupling_coefficient=one/dsqrt(bra_j_weight_product)
        two_bra_shell_j_3=bra_subshell_quantum_codes(3,orbital_index_3)-1
        two_ket_shell_j_3=ket_subshell_quantum_codes(3,orbital_index_3)-1
        recoupling_coefficient=recoupling_coefficient*dsqrt(dble(two_bra_shell_j_3+1))/dsqrt(dble((two_tensor_rank+1) &
            *(two_ket_shell_j_3+1)))
        coupling_allowed=0
        coupling_step_count=third_shell_position-second_shell_position
        if(coupling_step_count > 1) then
            call recouple_spectator_shell_chain(second_shell_position,third_shell_position,two_tensor_rank,evaluate_coefficient, &
                coupling_allowed,diagram_factor)
            if(coupling_allowed == 0)return
            recoupling_coefficient=diagram_factor*recoupling_coefficient
        end if
        coupling_allowed=0
        call recouple_later_active_shell(first_shell_position,third_shell_position,two_tensor_rank,evaluate_coefficient, &
            & coupling_allowed, &
            diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        coupling_allowed=0
        call recouple_pair_tensor(first_shell_position,second_shell_position,two_rank_1,two_rank_2,two_tensor_rank, &
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
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        return
    end subroutine recouple_three_subshells
end module three_subshell_recoupling
