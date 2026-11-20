! Reduce the one-to-two-electron tensor matrix element in quasispin.
module jj_reduced_cfp_two_electron_terms
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_two_electron_reduced_cfp
contains
    subroutine evaluate_two_electron_reduced_cfp(two_electron_term_index,one_electron_term_index,two_orbital_j,reduced_cfp)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, half, angular_tolerance
        use jj_subshell_term_decoding, only: decode_jj_subshell_term
        use quasispin_cg_rank_half, only: evaluate_quasispin_cg_rank_half
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: two_orbital_j, two_electron_term_index, one_electron_term_index
        real(real64), intent(out) :: reduced_cfp

        ! Local declarations
        integer      :: two_bra_q, two_bra_seniority, two_ket_q, two_ket_seniority, two_bra_term_j, two_ket_term_j, &
            bra_occupation, phase_exponent
        real(real64) :: normalization_squared, cg_coefficient, bra_q, operator_q_projection, ket_q, bra_q_projection, &
            ket_q_projection
        reduced_cfp=zero
        call decode_jj_subshell_term(two_electron_term_index,two_orbital_j,two_bra_q,two_bra_seniority,two_bra_term_j)
        call decode_jj_subshell_term(one_electron_term_index,two_orbital_j,two_ket_q,two_ket_seniority,two_ket_term_j)
        operator_q_projection=half
        bra_occupation=2
        bra_q=half*dble(two_bra_q)
        ket_q=half*dble(two_ket_q)
        bra_q_projection=-half*dble(half*(two_orbital_j+1)-bra_occupation)
        ket_q_projection=-half*dble(half*(two_orbital_j+1)-(bra_occupation-1))
        call evaluate_quasispin_cg_rank_half(ket_q,ket_q_projection,operator_q_projection,bra_q,bra_q_projection,cg_coefficient)
        if(dabs(cg_coefficient) < angular_tolerance) return
        normalization_squared=dble(bra_occupation*(two_bra_q+1)*(two_bra_term_j+1))
        reduced_cfp=dsqrt(normalization_squared)/cg_coefficient
        phase_exponent=-bra_occupation-1
        if((phase_exponent/2)*2 /= phase_exponent)reduced_cfp=-reduced_cfp
        return
    end subroutine evaluate_two_electron_reduced_cfp
end module jj_reduced_cfp_two_electron_terms
