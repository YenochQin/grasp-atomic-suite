! Build the completely reduced pair tensor from two reduced single-tensor elements.
module jj_reduced_w_from_cfp
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_reduced_w_from_cfp
contains
    subroutine evaluate_reduced_w_from_cfp(quasispin_rank,angular_rank,two_orbital_j,bra_term_index,ket_term_index, &
        & reduced_w_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, angular_tolerance
        use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j, first_opposite_parity_terms_low_j, &
            first_same_parity_terms_low_j
        use jj_subshell_term_decoding, only: decode_jj_subshell_term
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use jj_reduced_cfp, only: evaluate_jj_reduced_cfp
        use wigner_6j, only: wigner_6j_value
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: quasispin_rank, angular_rank, bra_term_index, ket_term_index
        integer,      intent(out) :: two_orbital_j
        real(real64), intent(out) :: reduced_w_element

        ! Local declarations
        integer :: two_quasispin_rank,two_angular_rank,intermediate_term_index,first_intermediate_term, &
            last_intermediate_term,two_bra_q,two_bra_seniority,two_bra_term_j,two_ket_q,two_ket_seniority, &
            two_ket_term_j, &
            two_intermediate_q,two_intermediate_seniority,two_intermediate_term_j
        real(real64) :: bra_reduced_cfp, ket_reduced_cfp, sum_over_terms, angular_six_j, quasispin_six_j
        reduced_w_element=zero
        if(first_same_parity_terms_low_j(bra_term_index) /= first_same_parity_terms_low_j(ket_term_index)) return
        sum_over_terms=zero
        call decode_jj_subshell_term(bra_term_index,two_orbital_j,two_bra_q,two_bra_seniority,two_bra_term_j)
        two_quasispin_rank=quasispin_rank*2
        two_angular_rank=angular_rank*2
        ! One single-electron tensor changes occupation by one, so the intermediate term
        ! belongs to the opposite occupation-parity set.
        first_intermediate_term=first_opposite_parity_terms_low_j(bra_term_index)
        last_intermediate_term=last_opposite_parity_terms_low_j(bra_term_index)
        call decode_jj_subshell_term(ket_term_index,two_orbital_j,two_ket_q,two_ket_seniority,two_ket_term_j)

        ! Recouple the two reduced single-tensor elements in angular and quasispin space.
        do intermediate_term_index=first_intermediate_term,last_intermediate_term
            call decode_jj_subshell_term(intermediate_term_index,two_orbital_j,two_intermediate_q,two_intermediate_seniority, &
                two_intermediate_term_j)
            if(allows_wigner_6j_triangles(two_orbital_j,two_orbital_j,two_angular_rank,two_ket_term_j,two_bra_term_j, &
                two_intermediate_term_j) /= 0) then
                if(allows_wigner_6j_triangles(1,1,two_quasispin_rank,two_ket_q,two_bra_q,two_intermediate_q) /= 0) then
                    call evaluate_jj_reduced_cfp(two_orbital_j,bra_term_index,two_bra_q,two_bra_term_j,intermediate_term_index, &
                        two_intermediate_q,two_intermediate_term_j,bra_reduced_cfp)
                    if(dabs(bra_reduced_cfp) > angular_tolerance) then
                        call evaluate_jj_reduced_cfp(two_orbital_j,intermediate_term_index,two_intermediate_q, &
                            & two_intermediate_term_j, &
                            ket_term_index,two_ket_q,two_ket_term_j,ket_reduced_cfp)
                        if(dabs(ket_reduced_cfp) > angular_tolerance) then
                            call wigner_6j_value(two_orbital_j,two_orbital_j,two_angular_rank,two_ket_term_j,two_bra_term_j, &
                                two_intermediate_term_j,0,angular_six_j)
                            call wigner_6j_value(1,1,two_quasispin_rank,two_ket_q,two_bra_q,two_intermediate_q,0,quasispin_six_j)
                            sum_over_terms=sum_over_terms+angular_six_j*quasispin_six_j*bra_reduced_cfp*ket_reduced_cfp
                        end if
                    end if
                end if
            end if
        end do
        reduced_w_element=sum_over_terms*dsqrt(dble((two_quasispin_rank+1)*(two_angular_rank+1)))
        if(mod(two_bra_q+two_bra_term_j+two_ket_q+two_ket_term_j+two_quasispin_rank+two_angular_rank, &
            4) /= 0)reduced_w_element=-reduced_w_element
        return
    end subroutine evaluate_reduced_w_from_cfp
end module jj_reduced_w_from_cfp
