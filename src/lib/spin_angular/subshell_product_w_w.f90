! Multiply two pair-tensor elements on distinct subshells.
module subshell_product_w_w
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_subshell_product_w_w
contains
    subroutine evaluate_subshell_product_w_w(first_w_rank,second_w_rank,q_projection_1,q_projection_2,q_projection_3, &
        & q_projection_4, &
        subshell_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, angular_tolerance
        use jj_tensor_coupling_tables, only: ket_shell_values_1, ket_shell_values_2, bra_shell_values_1, bra_shell_values_2, &
            ket_shell_codes_1, ket_shell_codes_2, bra_shell_codes_1, bra_shell_codes_2
        use subshell_tensor_diagnostics, only: report_subshell_tensor_error
        use jj_w_matrix_element, only: evaluate_jj_w
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)               :: first_w_rank, second_w_rank
        real(real64), intent(in)               :: q_projection_1, q_projection_2, q_projection_3, q_projection_4
        real(real64), intent(out)              :: subshell_matrix_element

        ! Local declarations
        integer :: occupation_change_1, occupation_change_2, occupation_change_3, occupation_change_4, &
            first_shell_occupation_change, last_shell_occupation_change
        real(real64) :: first_w_matrix_element, second_w_matrix_element
        subshell_matrix_element=zero
        if(bra_shell_codes_1(3) > 9) then
            if(bra_shell_codes_1(4) > 2) call report_subshell_tensor_error(35)
            if(ket_shell_codes_1(4) > 2) call report_subshell_tensor_error(35)
        end if
        if(bra_shell_codes_2(3) > 9) then
            if(bra_shell_codes_2(4) > 2) call report_subshell_tensor_error(35)
            if(ket_shell_codes_2(4) > 2) call report_subshell_tensor_error(35)
        end if
        occupation_change_1=q_projection_1+q_projection_1+tenth*q_projection_1
        occupation_change_2=q_projection_2+q_projection_2+tenth*q_projection_2
        first_shell_occupation_change=occupation_change_1+occupation_change_2
        if(bra_shell_codes_1(4) /= (ket_shell_codes_1(4)+first_shell_occupation_change))return
        occupation_change_3=q_projection_3+q_projection_3+tenth*q_projection_3
        occupation_change_4=q_projection_4+q_projection_4+tenth*q_projection_4
        last_shell_occupation_change=occupation_change_3+occupation_change_4
        if(bra_shell_codes_2(4) /= (ket_shell_codes_2(4)+last_shell_occupation_change))return
        call evaluate_jj_w(bra_shell_codes_1,bra_shell_values_1,ket_shell_codes_1,ket_shell_values_1,first_w_rank, &
            q_projection_1,q_projection_2,first_w_matrix_element)
        if(dabs(first_w_matrix_element) < angular_tolerance)return
        call evaluate_jj_w(bra_shell_codes_2,bra_shell_values_2,ket_shell_codes_2,ket_shell_values_2,second_w_rank, &
            q_projection_3,q_projection_4,second_w_matrix_element)
        if(dabs(second_w_matrix_element) < angular_tolerance)return
        subshell_matrix_element=first_w_matrix_element*second_w_matrix_element
        return
    end subroutine evaluate_subshell_product_w_w
end module subshell_product_w_w
