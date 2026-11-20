! Evaluate the 6j factor at the first active subshell.
module first_active_shell_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: diaga1
contains
    subroutine diaga1(first_shell_position,two_tensor_rank,evaluate_coefficient,coupling_allowed, &
        recoupling_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use csf_pair_angular_state,             only: peel_orbital_indices => jlist, bra_subshell_quantum_codes => jjq1, &
            ket_subshell_quantum_codes => jjq2, bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2
        use wigner_6j_triangle_selection, only: ixjtik
        use wigner_6j, only: sixj
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: first_shell_position, two_tensor_rank, evaluate_coefficient
        integer, intent(inout)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer :: orbital_index,two_bra_subshell_j,two_ket_subshell_j,coupling_position,two_spectator_j, &
            two_bra_coupled_j,two_ket_coupled_j,phase_sum,spectator_orbital
        real(real64) :: diagram_factor
        recoupling_coefficient = zero
        orbital_index=peel_orbital_indices(first_shell_position)
        two_bra_subshell_j=bra_subshell_quantum_codes(3,orbital_index)-1
        two_ket_subshell_j=ket_subshell_quantum_codes(3,orbital_index)-1
        if(first_shell_position <= 2) then
            spectator_orbital=peel_orbital_indices(1)
            if(first_shell_position == 1)spectator_orbital=peel_orbital_indices(2)
            two_spectator_j=bra_subshell_quantum_codes(3,spectator_orbital)-1
            two_bra_coupled_j=bra_intermediate_couplings(1)-1
            two_ket_coupled_j=ket_intermediate_couplings(1)-1
        else
            coupling_position=first_shell_position-2
            two_spectator_j=bra_intermediate_couplings(coupling_position)-1
            coupling_position=first_shell_position-1
            two_bra_coupled_j=bra_intermediate_couplings(coupling_position)-1
            two_ket_coupled_j=ket_intermediate_couplings(coupling_position)-1
        end if
        if(evaluate_coefficient /= 0) then
            call sixj(two_tensor_rank,two_ket_subshell_j,two_bra_subshell_j,two_spectator_j,two_bra_coupled_j, &
                two_ket_coupled_j,0,diagram_factor)
            diagram_factor=diagram_factor*dsqrt(dble((two_bra_subshell_j+1)*(two_ket_coupled_j+1)))
            phase_sum=two_spectator_j+two_bra_coupled_j+two_ket_subshell_j+two_tensor_rank
            if((phase_sum/4)*4 /= phase_sum)diagram_factor=-diagram_factor
            recoupling_coefficient=diagram_factor
            coupling_allowed=1

            ! Acting on the first member of the initial pair requires an interchange phase.
            if(first_shell_position /= 1)return
            phase_sum=two_bra_subshell_j+two_ket_subshell_j+2*two_spectator_j-two_bra_coupled_j-two_ket_coupled_j
            if((phase_sum/4)*4 /= phase_sum)recoupling_coefficient=-recoupling_coefficient
        else
            if(ixjtik(two_tensor_rank,two_ket_subshell_j,two_bra_subshell_j,two_spectator_j,two_bra_coupled_j, &
                two_ket_coupled_j) == 0)return
            coupling_allowed=1
        end if
        return
    end subroutine diaga1
end module first_active_shell_recoupling
