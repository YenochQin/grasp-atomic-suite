! Multiply the coupling-chain factors for two active subshells.
module two_subshell_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: reco2
contains
    subroutine reco2(first_shell_position,second_shell_position,two_tensor_rank,evaluate_coefficient, &
        coupling_allowed,recoupling_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: one
        use csf_pair_angular_state,             only: peel_orbital_indices => jlist, bra_subshell_quantum_codes => jjq1, &
            ket_subshell_quantum_codes => jjq2
        use first_active_shell_recoupling, only: diaga1
        use later_active_shell_recoupling, only: diaga2
        use spectator_shell_chain_recoupling, only: diaga3
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: first_shell_position,second_shell_position,two_tensor_rank,evaluate_coefficient
        integer, intent(out)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer      :: two_bra_shell_j_1,two_bra_shell_j_2,two_ket_shell_j_1,two_ket_shell_j_2,orbital_index_1, &
            orbital_index_2,coupling_step_count
        real(real64) :: first_bra_j_weight, bra_j_weight_product, diagram_factor
        coupling_allowed=0
        orbital_index_1=peel_orbital_indices(first_shell_position)
        orbital_index_2=peel_orbital_indices(second_shell_position)
        first_bra_j_weight=dble(bra_subshell_quantum_codes(3,orbital_index_1))
        bra_j_weight_product=dble(bra_subshell_quantum_codes(3,orbital_index_2))
        bra_j_weight_product=first_bra_j_weight*bra_j_weight_product
        recoupling_coefficient=one/dsqrt(bra_j_weight_product)
        if(evaluate_coefficient == 0) then
            coupling_allowed=0
        else if(two_tensor_rank /= 0) then
            coupling_allowed=0
        else
            coupling_allowed=1
            return
        end if
        two_bra_shell_j_1=bra_subshell_quantum_codes(3,orbital_index_1)-1
        two_bra_shell_j_2=bra_subshell_quantum_codes(3,orbital_index_2)-1
        two_ket_shell_j_1=ket_subshell_quantum_codes(3,orbital_index_1)-1
        two_ket_shell_j_2=ket_subshell_quantum_codes(3,orbital_index_2)-1
        call diaga2(first_shell_position,second_shell_position,two_tensor_rank,evaluate_coefficient,coupling_allowed, &
            diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient*dsqrt(dble(two_bra_shell_j_2+1))/ &
            dsqrt(dble((two_tensor_rank+1)*(two_ket_shell_j_2+1)))
        if(first_shell_position == 1.and.second_shell_position == 2)return
        coupling_allowed=0
        call diaga1(first_shell_position,two_tensor_rank,evaluate_coefficient,coupling_allowed,diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        coupling_step_count=second_shell_position-first_shell_position
        if(first_shell_position == 1)coupling_step_count=second_shell_position-1-first_shell_position
        if(coupling_step_count <= 1)return
        coupling_allowed=0
        call diaga3(first_shell_position,second_shell_position,two_tensor_rank,evaluate_coefficient,coupling_allowed, &
            diagram_factor)
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        return
    end subroutine reco2
end module two_subshell_recoupling
