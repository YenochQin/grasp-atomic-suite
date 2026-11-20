! Couple a pair tensor followed by a single creation or annihilation tensor.
module jj_wa_matrix_element
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_jj_wa
contains
    subroutine evaluate_jj_wa(bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,w_angular_rank, &
        total_angular_rank,q_projection_1,q_projection_2,q_projection_3,wa_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, two, angular_tolerance
        use subshell_tensor_diagnostics, only: report_subshell_tensor_error
        use jj_wa_matrix_element_high_j, only: evaluate_high_j_wa
        use quasispin_projection_selection, only: allows_quasispin_projection
        use jj_intermediate_terms_bra_parity, only: select_bra_parity_intermediate_terms
        use jj_subshell_term_decoding, only: decode_jj_subshell_term
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use quasispin_cg_rank_half, only: evaluate_quasispin_cg_rank_half
        use jj_reduced_cfp, only: evaluate_jj_reduced_cfp
        use jj_w_matrix_element, only: evaluate_jj_w
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)               :: w_angular_rank
        integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
        real(real64), intent(in)               :: total_angular_rank, q_projection_1, q_projection_2, q_projection_3
        real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
        real(real64), intent(out)              :: wa_matrix_element

        ! Local declarations
        integer :: two_w_rank, two_total_rank, twice_phase_exponent, occupation_change_3, &
            two_abs_intermediate_projection, intermediate_term_index, first_intermediate_term, &
            last_intermediate_term, two_intermediate_seniority
        integer :: intermediate_shell_codes(7)
        real(real64)               :: sum_over_terms, a_matrix_element_term, reduced_cfp, six_j_value, w_matrix_element
        real(real64) :: intermediate_shell_values(3)
        wa_matrix_element=zero
        if(ket_shell_codes(3) == 9) then
            if(max0(bra_shell_codes(4),ket_shell_codes(4)) < 3) then
                if(bra_shell_codes(1) < 300) call report_subshell_tensor_error(55)
                if(ket_shell_codes(1) < 300) call report_subshell_tensor_error(55)
                call evaluate_high_j_wa(w_angular_rank,total_angular_rank,q_projection_1,q_projection_2,q_projection_3, &
                    bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,wa_matrix_element)
                return
            end if
        else if(ket_shell_codes(3) > 9) then
            call evaluate_high_j_wa(w_angular_rank,total_angular_rank,q_projection_1,q_projection_2,q_projection_3, &
                bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,wa_matrix_element)
            return
        end if
        if(allows_quasispin_projection(ket_shell_codes(7),ket_shell_values(3),bra_shell_codes(7),bra_shell_values(3)) == 0)return
        sum_over_terms=zero
        two_w_rank=w_angular_rank*2
        occupation_change_3=q_projection_3*two
        two_total_rank=total_angular_rank+total_angular_rank+tenth*total_angular_rank
        if(select_bra_parity_intermediate_terms(bra_shell_codes,ket_shell_codes,two_total_rank,bra_shell_values,ket_shell_values, &
            intermediate_shell_codes,intermediate_shell_values,first_intermediate_term,last_intermediate_term, &
            occupation_change_3) == 0)return
        two_abs_intermediate_projection=two*dabs(intermediate_shell_values(3))+tenth

        ! Sum over admissible intermediate terms at the required occupation and projection.
        do intermediate_term_index=first_intermediate_term,last_intermediate_term
            call decode_jj_subshell_term(intermediate_term_index,intermediate_shell_codes(3),intermediate_shell_codes(7), &
                two_intermediate_seniority,intermediate_shell_codes(6))
            if(two_abs_intermediate_projection > intermediate_shell_codes(7)) cycle
            if(allows_wigner_6j_triangles(two_w_rank,bra_shell_codes(3),two_total_rank,ket_shell_codes(6),bra_shell_codes(6), &
                intermediate_shell_codes(6)) == 0) cycle
            intermediate_shell_codes(1)=intermediate_term_index
            intermediate_shell_values(2)=dble(intermediate_shell_codes(6))/two
            intermediate_shell_values(1)=dble(intermediate_shell_codes(7))/two
            call evaluate_quasispin_cg_rank_half(ket_shell_values(1),ket_shell_values(3),q_projection_3, &
                & intermediate_shell_values(1), &
                intermediate_shell_values(3),a_matrix_element_term)
            if(dabs(a_matrix_element_term) < angular_tolerance) cycle
            call evaluate_jj_reduced_cfp(ket_shell_codes(3),intermediate_shell_codes(1),intermediate_shell_codes(7), &
                intermediate_shell_codes(6),ket_shell_codes(1),ket_shell_codes(7),ket_shell_codes(6),reduced_cfp)
            if(dabs(reduced_cfp) < angular_tolerance) cycle
            a_matrix_element_term=a_matrix_element_term*reduced_cfp
            call evaluate_jj_w(bra_shell_codes,bra_shell_values,intermediate_shell_codes,intermediate_shell_values, &
                w_angular_rank,q_projection_1,q_projection_2,w_matrix_element)
            if(dabs(w_matrix_element) < angular_tolerance) cycle
            a_matrix_element_term=a_matrix_element_term*w_matrix_element
            call wigner_6j_value(two_w_rank,bra_shell_codes(3),two_total_rank,ket_shell_codes(6),bra_shell_codes(6), &
                intermediate_shell_codes(6),0,six_j_value)
            a_matrix_element_term=a_matrix_element_term*six_j_value/dsqrt(dble(intermediate_shell_codes(7)+1))
            sum_over_terms=sum_over_terms+a_matrix_element_term
        end do
        wa_matrix_element=sum_over_terms*dsqrt(dble(two_total_rank+1))
        twice_phase_exponent=two_total_rank+bra_shell_codes(6)+ket_shell_codes(6)+2
        if(((twice_phase_exponent/4)*4) /= twice_phase_exponent)wa_matrix_element=-wa_matrix_element
        return
    end subroutine evaluate_jj_wa
end module jj_wa_matrix_element
