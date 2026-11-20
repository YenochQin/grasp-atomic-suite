! Couple a pair tensor followed by a single creation or annihilation tensor.
module jj_wa_matrix_element
    implicit none (type, external)

    ! Module interface
    private
    public :: wap1
contains
    subroutine wap1(bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,w_angular_rank, &
        total_angular_rank,q_projection_1,q_projection_2,q_projection_3,wa_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, two, angular_tolerance => eps
        use subshell_tensor_diagnostics, only: mes
        use jj_wa_matrix_element_high_j, only: wap1jjg
        use quasispin_projection_selection, only: izas1
        use jj_intermediate_terms_bra_parity, only: itjj3
        use jj_subshell_term_decoding, only: rumtjj
        use wigner_6j_triangle_selection, only: ixjtik
        use quasispin_cg_rank_half, only: c0t5s
        use jj_reduced_cfp, only: rmeajj
        use jj_w_matrix_element, only: wj1
        use wigner_6j, only: sixj
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
                if(bra_shell_codes(1) < 300) call mes(55)
                if(ket_shell_codes(1) < 300) call mes(55)
                call wap1jjg(w_angular_rank,total_angular_rank,q_projection_1,q_projection_2,q_projection_3, &
                    bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,wa_matrix_element)
                return
            end if
        else if(ket_shell_codes(3) > 9) then
            call wap1jjg(w_angular_rank,total_angular_rank,q_projection_1,q_projection_2,q_projection_3, &
                bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,wa_matrix_element)
            return
        end if
        if(izas1(ket_shell_codes(7),ket_shell_values(3),bra_shell_codes(7),bra_shell_values(3)) == 0)return
        sum_over_terms=zero
        two_w_rank=w_angular_rank*2
        occupation_change_3=q_projection_3*two
        two_total_rank=total_angular_rank+total_angular_rank+tenth*total_angular_rank
        if(itjj3(bra_shell_codes,ket_shell_codes,two_total_rank,bra_shell_values,ket_shell_values, &
            intermediate_shell_codes,intermediate_shell_values,first_intermediate_term,last_intermediate_term, &
            occupation_change_3) == 0)return
        two_abs_intermediate_projection=two*dabs(intermediate_shell_values(3))+tenth

        ! Sum over admissible intermediate terms at the required occupation and projection.
        do intermediate_term_index=first_intermediate_term,last_intermediate_term
            call rumtjj(intermediate_term_index,intermediate_shell_codes(3),intermediate_shell_codes(7), &
                two_intermediate_seniority,intermediate_shell_codes(6))
            if(two_abs_intermediate_projection > intermediate_shell_codes(7)) cycle
            if(ixjtik(two_w_rank,bra_shell_codes(3),two_total_rank,ket_shell_codes(6),bra_shell_codes(6), &
                intermediate_shell_codes(6)) == 0) cycle
            intermediate_shell_codes(1)=intermediate_term_index
            intermediate_shell_values(2)=dble(intermediate_shell_codes(6))/two
            intermediate_shell_values(1)=dble(intermediate_shell_codes(7))/two
            call c0t5s(ket_shell_values(1),ket_shell_values(3),q_projection_3,intermediate_shell_values(1), &
                intermediate_shell_values(3),a_matrix_element_term)
            if(dabs(a_matrix_element_term) < angular_tolerance) cycle
            call rmeajj(ket_shell_codes(3),intermediate_shell_codes(1),intermediate_shell_codes(7), &
                intermediate_shell_codes(6),ket_shell_codes(1),ket_shell_codes(7),ket_shell_codes(6),reduced_cfp)
            if(dabs(reduced_cfp) < angular_tolerance) cycle
            a_matrix_element_term=a_matrix_element_term*reduced_cfp
            call wj1(bra_shell_codes,bra_shell_values,intermediate_shell_codes,intermediate_shell_values, &
                w_angular_rank,q_projection_1,q_projection_2,w_matrix_element)
            if(dabs(w_matrix_element) < angular_tolerance) cycle
            a_matrix_element_term=a_matrix_element_term*w_matrix_element
            call sixj(two_w_rank,bra_shell_codes(3),two_total_rank,ket_shell_codes(6),bra_shell_codes(6), &
                intermediate_shell_codes(6),0,six_j_value)
            a_matrix_element_term=a_matrix_element_term*six_j_value/dsqrt(dble(intermediate_shell_codes(7)+1))
            sum_over_terms=sum_over_terms+a_matrix_element_term
        end do
        wa_matrix_element=sum_over_terms*dsqrt(dble(two_total_rank+1))
        twice_phase_exponent=two_total_rank+bra_shell_codes(6)+ket_shell_codes(6)+2
        if(((twice_phase_exponent/4)*4) /= twice_phase_exponent)wa_matrix_element=-wa_matrix_element
        return
    end subroutine wap1
end module jj_wa_matrix_element
