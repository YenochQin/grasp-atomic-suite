! Form an off-diagonal one-electron tensor coefficient and its fermionic phase.
module one_electron_tensor_between_subshells
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_tensor_transfer_coefficient
contains
    subroutine calculate_tensor_transfer_coefficient(active_shell_count,tensor_rank,bra_shell_position,ket_shell_position, &
        angular_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, half, angular_tolerance
        use csf_pair_angular_state,             only: bra_occupations, peel_orbital_indices
        use orbital_csf_basis,           only: orbital_kappa
        use jj_tensor_coupling_tables, only: ket_shell_values_1, ket_shell_values_2, bra_shell_values_1, bra_shell_values_2, &
            ket_shell_codes_1, ket_shell_codes_2, bra_shell_codes_1, bra_shell_codes_2
        use one_electron_tensor_recoupling_selection, only: select_one_electron_tensor_recoupling
        use one_electron_two_subshell_recoupling, only: recouple_one_electron_two_subshells
        use quasispin_cg_rank_half, only: evaluate_quasispin_cg_rank_half
        use jj_reduced_cfp, only: evaluate_jj_reduced_cfp
        use jj_w_matrix_element, only: evaluate_jj_w
        use jj_active_subshell_states, only: prepare_jj_active_subshell_states
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: active_shell_count,tensor_rank,bra_shell_position,ket_shell_position
        real(real64), intent(out) :: angular_coefficient

        ! Local declarations
        integer      :: shell_position, first_shell_position, last_shell_position, coupling_allowed, orbital_index, phase_sum, &
            occupation_change_1, occupation_change_2, last_phase_shell, two_first_orbital_j, two_last_orbital_j
        real(real64) :: recoupling_factor,subshell_matrix_element,first_cg_factor,second_cg_factor,first_reduced_cfp, &
            second_reduced_cfp,q_projection_1,q_projection_2
        angular_coefficient=zero
        first_shell_position=min0(bra_shell_position,ket_shell_position)
        last_shell_position=max0(bra_shell_position,ket_shell_position)
        call select_one_electron_tensor_recoupling(active_shell_count,first_shell_position,last_shell_position,tensor_rank, &
            & coupling_allowed)
        if(coupling_allowed == 0)return
        orbital_index=peel_orbital_indices(first_shell_position)
        two_first_orbital_j=(iabs(orbital_kappa(orbital_index))*2)-1
        orbital_index=peel_orbital_indices(last_shell_position)
        two_last_orbital_j=(iabs(orbital_kappa(orbital_index))*2)-1
        call recouple_one_electron_two_subshells(active_shell_count,first_shell_position,last_shell_position,two_first_orbital_j, &
            two_last_orbital_j,tensor_rank,0,coupling_allowed,recoupling_factor)
        if(coupling_allowed == 0)return
        call prepare_jj_active_subshell_states(bra_shell_position,ket_shell_position,bra_shell_position,bra_shell_position,2)
        q_projection_1=half
        q_projection_2=-half
        occupation_change_1=q_projection_1+q_projection_1+tenth*q_projection_1
        if(bra_shell_codes_1(4) /= (ket_shell_codes_1(4)+occupation_change_1)) return
        occupation_change_2=q_projection_2+q_projection_2+tenth*q_projection_2
        if(bra_shell_codes_2(4) /= (ket_shell_codes_2(4)+occupation_change_2)) return
        call evaluate_quasispin_cg_rank_half(ket_shell_values_1(1),ket_shell_values_1(3),q_projection_1,bra_shell_values_1(1), &
            bra_shell_values_1(3),first_cg_factor)
        if(dabs(first_cg_factor) < angular_tolerance) return
        call evaluate_quasispin_cg_rank_half(ket_shell_values_2(1),ket_shell_values_2(3),q_projection_2,bra_shell_values_2(1), &
            bra_shell_values_2(3),second_cg_factor)
        if(dabs(second_cg_factor) < angular_tolerance) return
        call evaluate_jj_reduced_cfp(bra_shell_codes_1(3),bra_shell_codes_1(1),bra_shell_codes_1(7),bra_shell_codes_1(6), &
            ket_shell_codes_1(1),ket_shell_codes_1(7),ket_shell_codes_1(6),first_reduced_cfp)
        if(dabs(first_reduced_cfp) < angular_tolerance) return
        call evaluate_jj_reduced_cfp(bra_shell_codes_2(3),bra_shell_codes_2(1),bra_shell_codes_2(7),bra_shell_codes_2(6), &
            ket_shell_codes_2(1),ket_shell_codes_2(7),ket_shell_codes_2(6),second_reduced_cfp)
        if(dabs(second_reduced_cfp) < angular_tolerance) return
        subshell_matrix_element=first_reduced_cfp*second_reduced_cfp*first_cg_factor*second_cg_factor
        call recouple_one_electron_two_subshells(active_shell_count,first_shell_position,last_shell_position,two_first_orbital_j, &
            two_last_orbital_j,tensor_rank,1,coupling_allowed,recoupling_factor)
        angular_coefficient=subshell_matrix_element*recoupling_factor/dsqrt(dble((2*tensor_rank+1)* &
            (bra_shell_codes_1(7)+1)*(bra_shell_codes_2(7)+1)))
        last_phase_shell=last_shell_position-1
        phase_sum=0
        do  shell_position=first_shell_position,last_phase_shell
            orbital_index=peel_orbital_indices(shell_position)
            phase_sum=phase_sum+bra_occupations(orbital_index)
        end do
        phase_sum=phase_sum+1
        if(mod(phase_sum,2) /= 0)angular_coefficient=-angular_coefficient

        ! Reversing the transfer order contributes a separate angular-momentum phase.
        if(first_shell_position/=bra_shell_position) then
            phase_sum = bra_shell_codes_1(3)+bra_shell_codes_2(3)-2*tensor_rank+2
            if(mod(phase_sum,4) /= 0)angular_coefficient=-angular_coefficient
        end if
        angular_coefficient=-angular_coefficient*sqrt(dble(ket_shell_codes_1(3)+1))
        return
    end subroutine calculate_tensor_transfer_coefficient
end module one_electron_tensor_between_subshells
