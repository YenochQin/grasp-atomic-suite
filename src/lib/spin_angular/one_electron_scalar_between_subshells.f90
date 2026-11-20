! Form a scalar one-electron coefficient for transfer between subshells.
module one_electron_scalar_between_subshells
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_scalar_transfer_coefficient
contains
    subroutine calculate_scalar_transfer_coefficient(bra_csf_index,ket_csf_index,bra_shell_position,ket_shell_position, &
        & angular_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, half, angular_tolerance
        use csf_pair_angular_state, only: peel_orbital_indices, bra_occupations
        use jj_tensor_coupling_tables, only: ket_shell_values_1, ket_shell_values_2, bra_shell_values_1, bra_shell_values_2, &
            ket_shell_codes_1, ket_shell_codes_2, bra_shell_codes_1, bra_shell_codes_2
        use one_electron_scalar_recoupling_selection, only: select_one_electron_scalar_recoupling
        use jj_active_subshell_states, only: prepare_jj_active_subshell_states
        use two_subshell_recoupling, only: recouple_two_subshells
        use subshell_product_a_a, only: evaluate_subshell_product_a_a
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_position,ket_shell_position
        real(real64), intent(out) :: angular_coefficient

        ! Local declarations
        integer      :: coupling_allowed,first_shell_position,last_shell_position,intervening_electron_count, &
            last_phase_shell,shell_position,orbital_index,phase_sum
        real(real64) :: q_projection_1,q_projection_2,recoupling_factor,subshell_matrix_element
        angular_coefficient=zero
        if(bra_shell_position == ket_shell_position) return
        if(bra_shell_position < ket_shell_position) then
            first_shell_position=bra_shell_position
            last_shell_position=ket_shell_position
        else
            first_shell_position=ket_shell_position
            last_shell_position=bra_shell_position
        end if
        call select_one_electron_scalar_recoupling(-1,first_shell_position,last_shell_position,last_shell_position, &
            & last_shell_position,1, &
            coupling_allowed)
        if(coupling_allowed == 0)return
        q_projection_1=half
        q_projection_2=-half
        call prepare_jj_active_subshell_states(bra_shell_position,ket_shell_position,bra_shell_position,bra_shell_position,2)
        if(ket_shell_codes_1(3) /= ket_shell_codes_2(3)) return
        call recouple_two_subshells(first_shell_position,last_shell_position,ket_shell_codes_2(3),0,coupling_allowed, &
            & recoupling_factor)
        if(coupling_allowed == 0)return
        call evaluate_subshell_product_a_a(bra_shell_codes_1,bra_shell_codes_2,bra_shell_values_1,bra_shell_values_2, &
            & ket_shell_codes_1, &
            ket_shell_codes_2,ket_shell_values_1,ket_shell_values_2,q_projection_1,q_projection_2, &
            subshell_matrix_element)
        if(dabs(subshell_matrix_element) > angular_tolerance) then
            call recouple_two_subshells(first_shell_position,last_shell_position,ket_shell_codes_2(3),1,coupling_allowed, &
                recoupling_factor)
            angular_coefficient=subshell_matrix_element*recoupling_factor*dsqrt(dble(ket_shell_codes_1(3)+1))

            ! The transfer sign depends on the bra occupations between the ordered subshells.
            intervening_electron_count=0
            last_phase_shell=last_shell_position-1
            do shell_position=first_shell_position,last_phase_shell
                orbital_index=peel_orbital_indices(shell_position)
                intervening_electron_count=bra_occupations(orbital_index)+intervening_electron_count
            end do
            if((intervening_electron_count/2)*2 == intervening_electron_count) angular_coefficient=-angular_coefficient
            angular_coefficient=-angular_coefficient
            phase_sum=bra_shell_codes_1(5)*bra_shell_codes_1(4)+bra_shell_codes_2(5)*bra_shell_codes_2(4)- &
                ket_shell_codes_1(5)*ket_shell_codes_1(4)-ket_shell_codes_2(5)*ket_shell_codes_2(4)
            if((phase_sum/4)*4 /= phase_sum)angular_coefficient=-angular_coefficient
        end if
        return
    end subroutine calculate_scalar_transfer_coefficient
end module one_electron_scalar_between_subshells
