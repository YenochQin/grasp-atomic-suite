! Combine the WA product on one subshell with a single tensor on another.
module subshell_product_wa_a
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_subshell_product_wa_a
contains
    subroutine evaluate_subshell_product_wa_a(bra_shell_codes_1,bra_shell_codes_2,bra_shell_values_1,bra_shell_values_2, &
        ket_shell_codes_1,ket_shell_codes_2,ket_shell_values_1,ket_shell_values_2,w_angular_rank,  &
        q_projection_1,q_projection_2,q_projection_3,q_projection_4,subshell_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, half, angular_tolerance
        use subshell_tensor_diagnostics, only: report_subshell_tensor_error
        use quasispin_cg_rank_half, only: evaluate_quasispin_cg_rank_half
        use jj_reduced_cfp, only: evaluate_jj_reduced_cfp
        use jj_wa_matrix_element, only: evaluate_jj_wa
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)               :: w_angular_rank
        integer,      intent(in) :: bra_shell_codes_1(7), bra_shell_codes_2(7), ket_shell_codes_1(7), &
            ket_shell_codes_2(7)
        real(real64), intent(in)               :: q_projection_1, q_projection_2, q_projection_3, q_projection_4
        real(real64), intent(in) :: bra_shell_values_1(3), bra_shell_values_2(3), ket_shell_values_1(3), &
            ket_shell_values_2(3)
        real(real64), intent(out)              :: subshell_matrix_element

        ! Local declarations
        integer :: occupation_change_1, occupation_change_2, occupation_change_3, occupation_change_4, &
            first_shell_occupation_change, two_w_rank
        real(real64) :: cg_factor, reduced_cfp, three_operator_rank, wa_matrix_element
        subshell_matrix_element=zero
        if(bra_shell_codes_1(3) > 9) then
            if(bra_shell_codes_1(4) > 2) call report_subshell_tensor_error(34)
            if(ket_shell_codes_1(4) > 2) call report_subshell_tensor_error(34)
        end if
        if(bra_shell_codes_2(3) > 9) then
            if(bra_shell_codes_2(4) > 2) call report_subshell_tensor_error(34)
            if(ket_shell_codes_2(4) > 2) call report_subshell_tensor_error(34)
        end if
        occupation_change_1=q_projection_1+q_projection_1+tenth*q_projection_1
        occupation_change_2=q_projection_2+q_projection_2+tenth*q_projection_2
        occupation_change_3=q_projection_3+q_projection_3+tenth*q_projection_3
        first_shell_occupation_change=occupation_change_1+occupation_change_2+occupation_change_3
        if(bra_shell_codes_1(4) /= (ket_shell_codes_1(4)+first_shell_occupation_change))return
        occupation_change_4=q_projection_4+q_projection_4+tenth*q_projection_4
        if(bra_shell_codes_2(4) /= (ket_shell_codes_2(4)+occupation_change_4))return
        two_w_rank=w_angular_rank*2
        call evaluate_quasispin_cg_rank_half(ket_shell_values_2(1),ket_shell_values_2(3),q_projection_4,bra_shell_values_2(1), &
            bra_shell_values_2(3),cg_factor)
        if(abs(cg_factor) < angular_tolerance)return
        call evaluate_jj_reduced_cfp(bra_shell_codes_2(3),bra_shell_codes_2(1),bra_shell_codes_2(7),bra_shell_codes_2(6), &
            ket_shell_codes_2(1),ket_shell_codes_2(7),ket_shell_codes_2(6),reduced_cfp)
        if(abs(reduced_cfp) < angular_tolerance)return
        three_operator_rank=half*dble(bra_shell_codes_2(3))
        call evaluate_jj_wa(bra_shell_codes_1,bra_shell_values_1,ket_shell_codes_1,ket_shell_values_1,w_angular_rank, &
            three_operator_rank,q_projection_1,q_projection_2,q_projection_3,wa_matrix_element)
        if(abs(wa_matrix_element) < angular_tolerance)return
        subshell_matrix_element=-cg_factor*wa_matrix_element*reduced_cfp/sqrt(dble(bra_shell_codes_2(7)+1))
        return
    end subroutine evaluate_subshell_product_wa_a
end module subshell_product_wa_a
