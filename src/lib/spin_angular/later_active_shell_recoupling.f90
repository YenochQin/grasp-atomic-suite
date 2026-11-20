! Evaluate the 6j factor where a tensor acts on a later subshell.
module later_active_shell_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: recouple_later_active_shell
contains
    subroutine recouple_later_active_shell(first_shell_position,second_shell_position,two_tensor_rank,evaluate_coefficient, &
        coupling_allowed,recoupling_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use csf_pair_angular_state,             only: peel_orbital_indices, bra_subshell_quantum_codes, &
            ket_subshell_quantum_codes, bra_intermediate_couplings, ket_intermediate_couplings
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: first_shell_position, second_shell_position, two_tensor_rank, evaluate_coefficient
        integer, intent(inout)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer :: orbital_index_1,orbital_index_2,two_bra_shell_j_1,two_bra_shell_j_2,two_ket_shell_j_1, &
            two_ket_shell_j_2,coupled_position,preceding_coupling_position,two_bra_coupled_j,two_bra_prefix_j, &
            two_ket_prefix_j,phase_sum
        real(real64) :: six_j_value
        recoupling_coefficient = zero
        orbital_index_1=peel_orbital_indices(first_shell_position)
        orbital_index_2=peel_orbital_indices(second_shell_position)
        two_bra_shell_j_1=bra_subshell_quantum_codes(3,orbital_index_1)-1
        two_bra_shell_j_2=bra_subshell_quantum_codes(3,orbital_index_2)-1
        two_ket_shell_j_1=ket_subshell_quantum_codes(3,orbital_index_1)-1
        two_ket_shell_j_2=ket_subshell_quantum_codes(3,orbital_index_2)-1
        if(first_shell_position == 1.and.second_shell_position == 2) then
            two_bra_prefix_j=two_bra_shell_j_1
            two_ket_prefix_j=two_ket_shell_j_1
            two_bra_coupled_j=bra_intermediate_couplings(1)-1
        else
            coupled_position=second_shell_position-1
            two_bra_coupled_j=bra_intermediate_couplings(coupled_position)-1
            preceding_coupling_position=second_shell_position-2
            two_bra_prefix_j=bra_intermediate_couplings(preceding_coupling_position)-1
            two_ket_prefix_j=ket_intermediate_couplings(preceding_coupling_position)-1
        end if
        if(evaluate_coefficient /= 0) then
            call wigner_6j_value(two_tensor_rank,two_ket_shell_j_2,two_bra_shell_j_2,two_bra_coupled_j,two_bra_prefix_j, &
                two_ket_prefix_j,0,six_j_value)
            recoupling_coefficient=six_j_value*dsqrt(dble((two_ket_shell_j_2+1)*(two_bra_prefix_j+1)))
            phase_sum=two_bra_coupled_j+two_ket_prefix_j+two_bra_shell_j_2+two_tensor_rank
            if((phase_sum/4)*4 /= phase_sum)recoupling_coefficient=-recoupling_coefficient
            coupling_allowed=1
        else
            if(allows_wigner_6j_triangles(two_tensor_rank,two_ket_shell_j_2,two_bra_shell_j_2,two_bra_coupled_j,two_bra_prefix_j, &
                two_ket_prefix_j) == 0)return
            coupling_allowed=1
        end if
        return
    end subroutine recouple_later_active_shell
end module later_active_shell_recoupling
