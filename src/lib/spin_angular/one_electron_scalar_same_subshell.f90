! Form a rank-zero one-electron coefficient on one subshell.
module one_electron_scalar_same_subshell
    implicit none (type, external)

    ! Module interface
    private
    public :: calculate_scalar_same_subshell_coefficient
contains
    subroutine calculate_scalar_same_subshell_coefficient(active_shell_count,bra_csf_index,ket_csf_index,bra_shell_position, &
        & ket_shell_position, &
        angular_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, half, one, angular_tolerance
        use jj_tensor_coupling_tables, only: ket_shell_values_1, bra_shell_values_1, ket_shell_codes_1, bra_shell_codes_1
        use one_electron_scalar_recoupling_selection, only: select_one_electron_scalar_recoupling
        use jj_active_subshell_states, only: prepare_jj_active_subshell_states
        use jj_w_matrix_element, only: evaluate_jj_w
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: active_shell_count,bra_csf_index,ket_csf_index,bra_shell_position, &
            ket_shell_position
        real(real64), intent(out) :: angular_coefficient

        ! Local declarations
        integer      :: coupling_allowed
        real(real64) :: w_matrix_element,q_projection_1,q_projection_2,recoupling_factor
        angular_coefficient=zero
        if(bra_shell_position == ket_shell_position) then
            if(bra_csf_index /= ket_csf_index) then
                call select_one_electron_scalar_recoupling(active_shell_count,bra_shell_position,bra_shell_position, &
                    & bra_shell_position, &
                    bra_shell_position,0,coupling_allowed)
                if(coupling_allowed == 0)return
            end if
            call prepare_jj_active_subshell_states(bra_shell_position,bra_shell_position,bra_shell_position,bra_shell_position,1)
            q_projection_1=half
            q_projection_2=-half
            call evaluate_jj_w(bra_shell_codes_1,bra_shell_values_1,ket_shell_codes_1,ket_shell_values_1,0,q_projection_1, &
                q_projection_2,w_matrix_element)
            if(dabs(w_matrix_element) > angular_tolerance) then
                recoupling_factor=one/dsqrt(dble(bra_shell_codes_1(6)+1))
                angular_coefficient=w_matrix_element*recoupling_factor*dsqrt(dble(ket_shell_codes_1(3)+1))
                angular_coefficient=-angular_coefficient
            end if
        end if
        return
    end subroutine calculate_scalar_same_subshell_coefficient
end module one_electron_scalar_same_subshell
