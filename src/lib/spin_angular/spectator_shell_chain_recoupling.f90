! Multiply 6j factors along spectator subshells between active positions.
module spectator_shell_chain_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: recouple_spectator_shell_chain
contains
    subroutine recouple_spectator_shell_chain(first_shell_position,second_shell_position,two_tensor_rank,evaluate_coefficient, &
        coupling_allowed,recoupling_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, one
        use csf_pair_angular_state,             only: peel_orbital_indices, bra_subshell_quantum_codes, &
            bra_intermediate_couplings, ket_intermediate_couplings
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: first_shell_position, second_shell_position, two_tensor_rank, evaluate_coefficient
        integer, intent(inout)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer :: shell_position,orbital_index,two_spectator_j,coupled_position,preceding_coupling_position, &
            two_bra_input_j,two_bra_coupled_j,two_ket_input_j,two_ket_coupled_j,phase_sum
        real(real64) :: diagram_product, diagram_factor
        recoupling_coefficient = zero

        ! The first two shells share the first coupling node, so skip that node when starting at shell 1.
        shell_position=first_shell_position+1
        if(first_shell_position == 1)shell_position=shell_position+1
        if(evaluate_coefficient == 0) then
            if(shell_position >= second_shell_position) then
                recoupling_coefficient=one
                coupling_allowed=1
                return
            end if
            do
                orbital_index=peel_orbital_indices(shell_position)
                two_spectator_j=bra_subshell_quantum_codes(3,orbital_index)-1
                preceding_coupling_position=shell_position-2
                two_bra_input_j=bra_intermediate_couplings(preceding_coupling_position)-1
                two_ket_input_j=ket_intermediate_couplings(preceding_coupling_position)-1
                coupled_position=shell_position-1
                two_bra_coupled_j=bra_intermediate_couplings(coupled_position)-1
                two_ket_coupled_j=ket_intermediate_couplings(coupled_position)-1
                if(allows_wigner_6j_triangles(two_tensor_rank,two_ket_input_j,two_bra_input_j,two_spectator_j,two_bra_coupled_j, &
                    two_ket_coupled_j) == 0)return
                shell_position=shell_position+1
                if(shell_position == second_shell_position) exit
            end do
            recoupling_coefficient=one
            coupling_allowed=1
            return
        end if
        diagram_product=one
        do while (shell_position < second_shell_position)
            orbital_index=peel_orbital_indices(shell_position)
            two_spectator_j=bra_subshell_quantum_codes(3,orbital_index)-1
            preceding_coupling_position=shell_position-2
            two_bra_input_j=bra_intermediate_couplings(preceding_coupling_position)-1
            two_ket_input_j=ket_intermediate_couplings(preceding_coupling_position)-1
            coupled_position=shell_position-1
            two_bra_coupled_j=bra_intermediate_couplings(coupled_position)-1
            two_ket_coupled_j=ket_intermediate_couplings(coupled_position)-1
            call wigner_6j_value(two_tensor_rank,two_ket_input_j,two_bra_input_j,two_spectator_j,two_bra_coupled_j, &
                two_ket_coupled_j,0,diagram_factor)
            diagram_factor=diagram_factor*sqrt(dble((two_bra_input_j+1)*(two_ket_coupled_j+1)))
            phase_sum=two_tensor_rank+two_spectator_j+two_bra_input_j+two_ket_coupled_j
            if((phase_sum/4)*4 /= phase_sum)diagram_factor=-diagram_factor
            diagram_product=diagram_product*diagram_factor
            shell_position=shell_position+1
        end do
        recoupling_coefficient=diagram_product
        coupling_allowed=1
        return
    end subroutine recouple_spectator_shell_chain
end module spectator_shell_chain_recoupling
