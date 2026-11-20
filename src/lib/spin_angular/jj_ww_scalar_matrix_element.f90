! Couple two equal-rank pair tensors to a scalar on one subshell.
module jj_ww_scalar_matrix_element
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_jj_scalar_ww
contains
    subroutine evaluate_jj_scalar_ww(bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,w_angular_rank, &
        q_projection_1,q_projection_2,q_projection_3,q_projection_4,ww_scalar_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, two, angular_tolerance
        use quasispin_projection_selection, only: allows_quasispin_projection
        use jj_intermediate_terms_same_parity, only: select_same_parity_intermediate_terms
        use jj_subshell_term_decoding, only: decode_jj_subshell_term
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use jj_w_matrix_element, only: evaluate_jj_w
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)               :: w_angular_rank
        integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
        real(real64), intent(in)               :: q_projection_1, q_projection_2, q_projection_3, q_projection_4
        real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
        real(real64), intent(out)              :: ww_scalar_matrix_element

        ! Local declarations
        integer :: two_w_rank,two_total_rank,intermediate_occupation_change,two_abs_intermediate_projection, &
            occupation_change_3,occupation_change_4,twice_phase_exponent,partial_phase_sum,intermediate_term_index, &
            first_intermediate_term,last_intermediate_term,two_intermediate_seniority
        integer      :: intermediate_shell_codes(7)
        real(real64)               :: sum_over_terms, bra_w_matrix_element, ket_w_matrix_element
        real(real64) :: intermediate_shell_values(3)
        ww_scalar_matrix_element=zero
        if(ket_shell_codes(6) /= bra_shell_codes(6))return
        if(allows_quasispin_projection(ket_shell_codes(7),ket_shell_values(3),bra_shell_codes(7),bra_shell_values(3)) == 0)return
        sum_over_terms=zero
        two_w_rank=w_angular_rank*2
        occupation_change_3=q_projection_3*two+q_projection_3*tenth
        occupation_change_4=q_projection_4*two+q_projection_4*tenth
        intermediate_occupation_change=occupation_change_3+occupation_change_4
        if(select_same_parity_intermediate_terms(bra_shell_codes,ket_shell_codes,0,bra_shell_values,ket_shell_values, &
            & intermediate_shell_codes, &
            intermediate_shell_values,two_total_rank,first_intermediate_term,last_intermediate_term, &
            intermediate_occupation_change) == 0)return
        partial_phase_sum=two_w_rank-bra_shell_codes(6)
        two_abs_intermediate_projection=two*dabs(intermediate_shell_values(3))+tenth

        ! Sum over admissible intermediate terms at the required occupation and projection.
        do intermediate_term_index=first_intermediate_term,last_intermediate_term
            call decode_jj_subshell_term(intermediate_term_index,intermediate_shell_codes(3),intermediate_shell_codes(7), &
                two_intermediate_seniority,intermediate_shell_codes(6))
            if(two_abs_intermediate_projection > intermediate_shell_codes(7)) cycle
            if(allows_wigner_6j_triangles(two_w_rank,two_w_rank,0,ket_shell_codes(6),bra_shell_codes(6), &
                intermediate_shell_codes(6)) == 0) cycle
            intermediate_shell_codes(1)=intermediate_term_index
            intermediate_shell_values(2)=dble(intermediate_shell_codes(6))/two
            intermediate_shell_values(1)=dble(intermediate_shell_codes(7))/two
            call evaluate_jj_w(bra_shell_codes,bra_shell_values,intermediate_shell_codes,intermediate_shell_values, &
                w_angular_rank,q_projection_1,q_projection_2,bra_w_matrix_element)
            if(dabs(bra_w_matrix_element) < angular_tolerance) cycle
            call evaluate_jj_w(intermediate_shell_codes,intermediate_shell_values,ket_shell_codes,ket_shell_values, &
                w_angular_rank,q_projection_3,q_projection_4,ket_w_matrix_element)
            if(dabs(ket_w_matrix_element) < angular_tolerance) cycle
            bra_w_matrix_element=bra_w_matrix_element*ket_w_matrix_element
            twice_phase_exponent=partial_phase_sum+intermediate_shell_codes(6)
            if(((twice_phase_exponent/4)*4) /= twice_phase_exponent)bra_w_matrix_element=-bra_w_matrix_element
            sum_over_terms=sum_over_terms+bra_w_matrix_element
        end do
        ww_scalar_matrix_element=sum_over_terms/dsqrt(dble(two_w_rank+1)*(bra_shell_codes(6)+1))
        return
    end subroutine evaluate_jj_scalar_ww
end module jj_ww_scalar_matrix_element
