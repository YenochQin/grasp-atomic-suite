! Evaluate the high-j AW product by summing over intermediate terms.
module jj_aw_matrix_element_high_j
    implicit none (type, external)

    ! Module interface
    private
    public :: awp1jjg
contains
    subroutine awp1jjg(w_angular_rank,total_angular_rank,q_projection_1,q_projection_2,q_projection_3, &
        bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,aw_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, two, tenth, angular_tolerance => eps
        use jj_intermediate_terms_ket_parity, only: itjj2
        use wigner_6j_triangle_selection, only: ixjtik
        use quasispin_projection_selection, only: izas1
        use jj_subshell_term_decoding, only: rumtjj
        use jj_a_matrix_element_few_electrons, only: a1jj
        use jj_w_matrix_element_high_j, only: w1jjg
        use wigner_6j, only: sixj
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)               :: w_angular_rank
        integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
        real(real64), intent(in)               :: total_angular_rank, q_projection_1, q_projection_2, q_projection_3
        real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
        real(real64), intent(out)              :: aw_matrix_element

        ! Local declarations
        integer :: two_w_rank, two_total_rank, twice_phase_exponent, intermediate_occupation_change, &
            occupation_change_2, occupation_change_3, two_abs_intermediate_projection, intermediate_term_index, &
            first_intermediate_term, last_intermediate_term, two_intermediate_seniority
        integer :: intermediate_shell_codes(7)
        real(real64)               :: sum_over_terms, a_matrix_element_term, six_j_value, w_matrix_element
        real(real64) :: intermediate_shell_values(3)
        aw_matrix_element=zero
        if(izas1(ket_shell_codes(7),ket_shell_values(3),bra_shell_codes(7),bra_shell_values(3)) == 0)return
        sum_over_terms=zero
        two_w_rank=w_angular_rank*2
        two_total_rank=total_angular_rank+total_angular_rank+tenth*total_angular_rank
        occupation_change_2=q_projection_2*two+q_projection_2*tenth
        occupation_change_3=q_projection_3*two+q_projection_3*tenth
        intermediate_occupation_change=occupation_change_2+occupation_change_2
        if(ket_shell_codes(3) > 37) return
        if(itjj2(bra_shell_codes,ket_shell_codes,two_total_rank,bra_shell_values,ket_shell_values, &
            intermediate_shell_codes,intermediate_shell_values,first_intermediate_term,last_intermediate_term, &
            intermediate_occupation_change) == 0)return
        two_abs_intermediate_projection=two*dabs(intermediate_shell_values(3))+tenth

        ! Sum over admissible intermediate terms at the required occupation and projection.
        do intermediate_term_index=first_intermediate_term,last_intermediate_term
            call rumtjj(intermediate_term_index,intermediate_shell_codes(3),intermediate_shell_codes(7), &
                two_intermediate_seniority,intermediate_shell_codes(6))
            if(two_abs_intermediate_projection <= intermediate_shell_codes(7)) then
                if(ixjtik(bra_shell_codes(3),two_w_rank,two_total_rank,ket_shell_codes(6),bra_shell_codes(6), &
                    intermediate_shell_codes(6)) /= 0) then
                    intermediate_shell_codes(1)=intermediate_term_index
                    intermediate_shell_values(2)=dble(intermediate_shell_codes(6))/two
                    intermediate_shell_values(1)=dble(intermediate_shell_codes(7))/two
                    call a1jj(bra_shell_codes,bra_shell_values,intermediate_shell_codes,intermediate_shell_values, &
                        q_projection_1,a_matrix_element_term)
                    if(dabs(a_matrix_element_term) > angular_tolerance) then
                        call w1jjg(w_angular_rank,q_projection_2,q_projection_3,intermediate_shell_codes, &
                            intermediate_shell_values,ket_shell_codes,ket_shell_values,w_matrix_element)
                        if(dabs(w_matrix_element) > angular_tolerance) then
                            a_matrix_element_term=a_matrix_element_term*w_matrix_element
                            call sixj(bra_shell_codes(3),two_w_rank,two_total_rank,ket_shell_codes(6),bra_shell_codes(6), &
                                intermediate_shell_codes(6),0,six_j_value)
                            sum_over_terms=sum_over_terms+a_matrix_element_term*six_j_value
                        end if
                    end if
                end if
            end if
        end do
        aw_matrix_element=sum_over_terms*dsqrt(dble((two_total_rank+1)))
        twice_phase_exponent=two_total_rank+bra_shell_codes(6)+ket_shell_codes(6)
        if(((twice_phase_exponent/4)*4) /= twice_phase_exponent)aw_matrix_element=-aw_matrix_element
        return
    end subroutine awp1jjg
end module jj_aw_matrix_element_high_j
