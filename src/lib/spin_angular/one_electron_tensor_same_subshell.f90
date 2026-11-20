! Form a one-electron tensor coefficient for unchanged occupations at one subshell.
module one_electron_tensor_same_subshell
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_tensor_same_subshell_coefficient
contains
    subroutine calculate_tensor_same_subshell_coefficient(active_shell_count,tensor_rank,bra_csf_index,ket_csf_index, &
        & bra_shell_position, &
        ket_shell_position,angular_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, half, angular_tolerance
        use jj_tensor_coupling_tables, only: ket_shell_values_1, bra_shell_values_1, ket_shell_codes_1, bra_shell_codes_1
        use one_electron_tensor_recoupling_selection, only: select_one_electron_tensor_recoupling
        use one_electron_same_subshell_recoupling, only: recouple_one_electron_same_subshell
        use one_electron_two_subshell_recoupling, only: recouple_one_electron_two_subshells
        use jj_w_matrix_element, only: evaluate_jj_w
        use jj_active_subshell_states, only: prepare_jj_active_subshell_states
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: active_shell_count,tensor_rank,bra_csf_index,ket_csf_index,bra_shell_position, &
            ket_shell_position
        real(real64), intent(out) :: angular_coefficient

        ! Local declarations
        integer      :: coupling_allowed
        real(real64) :: recoupling_factor,w_matrix_element,q_projection_1,q_projection_2
        angular_coefficient=zero
        if(bra_shell_position == ket_shell_position) then
            if(bra_csf_index /= ket_csf_index) then
                call select_one_electron_tensor_recoupling(active_shell_count,bra_shell_position,bra_shell_position,tensor_rank, &
                    & coupling_allowed)
                if(coupling_allowed == 0)return
            end if
            call recouple_one_electron_same_subshell(active_shell_count,bra_shell_position,tensor_rank,0,coupling_allowed, &
                & recoupling_factor)
            if(coupling_allowed == 0)return
            call prepare_jj_active_subshell_states(bra_shell_position,bra_shell_position,bra_shell_position,bra_shell_position,1)
            q_projection_1=half
            q_projection_2=-half
            call evaluate_jj_w(bra_shell_codes_1,bra_shell_values_1,ket_shell_codes_1,ket_shell_values_1,tensor_rank, &
                q_projection_1,q_projection_2,w_matrix_element)
            if(dabs(w_matrix_element) > angular_tolerance) then
                call recouple_one_electron_same_subshell(active_shell_count,bra_shell_position,tensor_rank,1,coupling_allowed, &
                    & recoupling_factor)
                angular_coefficient=w_matrix_element*recoupling_factor*dsqrt(dble(ket_shell_codes_1(3)+1))/dsqrt(dble(2* &
                    tensor_rank+1))
                angular_coefficient=-angular_coefficient
            end if
        end if
        return
    end subroutine calculate_tensor_same_subshell_coefficient
end module one_electron_tensor_same_subshell
