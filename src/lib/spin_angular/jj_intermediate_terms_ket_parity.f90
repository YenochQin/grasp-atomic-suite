! Select intermediate terms with the ket occupation parity.
module jj_intermediate_terms_ket_parity
    implicit none (type, external)

    ! Module interface
    private
    public :: select_ket_parity_intermediate_terms
contains
    integer function select_ket_parity_intermediate_terms(bra_shell_codes,ket_shell_codes,two_tensor_rank,bra_shell_values, &
        & ket_shell_values, &
        intermediate_shell_codes,intermediate_shell_values,first_intermediate_term,last_intermediate_term, &
        occupation_change) result(coupling_allowed)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: half
        use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j, last_same_parity_terms_low_j, &
            first_opposite_parity_terms_low_j, first_same_parity_terms_low_j
        use jj_reduced_coefficients_9, only: last_opposite_parity_terms_two_j_9, last_same_parity_terms_two_j_9, &
            first_opposite_parity_terms_two_j_9, first_same_parity_terms_two_j_9
        use jj_reduced_coefficients_11, only: last_opposite_parity_terms_high_j, last_same_parity_terms_high_j, &
            first_opposite_parity_terms_high_j, first_same_parity_terms_high_j
        use subshell_tensor_diagnostics, only: report_subshell_tensor_error
        use angular_triangle_and_parity, only: satisfies_angular_triangle_and_parity
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)                :: two_tensor_rank, occupation_change
        integer,      intent(out)               :: first_intermediate_term, last_intermediate_term
        integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
        integer,      intent(out) :: intermediate_shell_codes(7)
        real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
        real(real64), intent(out) :: intermediate_shell_values(3)

        ! Local declarations
        integer :: bra_term_lookup_index, ket_term_lookup_index, first_bra_range_term, last_bra_range_term
        coupling_allowed=0
        if(ket_shell_codes(3) > 37) return
        if(satisfies_angular_triangle_and_parity(ket_shell_codes(6),bra_shell_codes(6),two_tensor_rank) == 0)return
        bra_term_lookup_index=bra_shell_codes(1)
        ket_term_lookup_index=ket_shell_codes(1)
        if(ket_shell_codes(3) < 9) then
            first_bra_range_term=first_opposite_parity_terms_low_j(bra_term_lookup_index)
            first_intermediate_term=first_same_parity_terms_low_j(ket_term_lookup_index)

            ! Only a common term range can contribute to the intermediate sum.
            if(first_bra_range_term/=first_intermediate_term)return
            last_bra_range_term=last_opposite_parity_terms_low_j(bra_term_lookup_index)
            last_intermediate_term=last_same_parity_terms_low_j(ket_term_lookup_index)
        else if(ket_shell_codes(3) == 9) then
            if(bra_term_lookup_index > 300) then
                if(ket_term_lookup_index < 300) call report_subshell_tensor_error(52)
                if(ket_shell_codes(4) > 2) call report_subshell_tensor_error(12)
                if(bra_shell_codes(4) > 2) call report_subshell_tensor_error(12)
                bra_term_lookup_index=bra_term_lookup_index-300
                ket_term_lookup_index=ket_term_lookup_index-300
                first_bra_range_term=first_opposite_parity_terms_two_j_9(bra_term_lookup_index)
                first_intermediate_term=first_same_parity_terms_two_j_9(ket_term_lookup_index)
                if(first_bra_range_term /= first_intermediate_term)return
                last_bra_range_term=last_opposite_parity_terms_two_j_9(bra_term_lookup_index)
                last_intermediate_term=last_same_parity_terms_two_j_9(ket_term_lookup_index)
            else
                print*, "ERROR in select_ket_parity_intermediate_terms"
                stop
            end if
        else
            if(ket_shell_codes(4) > 2) call report_subshell_tensor_error(12)
            if(bra_shell_codes(4) > 2) call report_subshell_tensor_error(12)
            first_bra_range_term=first_opposite_parity_terms_high_j(bra_term_lookup_index)
            first_intermediate_term=first_same_parity_terms_high_j(ket_term_lookup_index)
            if(first_bra_range_term /= first_intermediate_term)return
            last_bra_range_term=last_opposite_parity_terms_high_j(bra_term_lookup_index)
            last_intermediate_term=last_same_parity_terms_high_j(ket_term_lookup_index)
        end if
        if(last_bra_range_term /= last_intermediate_term)return
        coupling_allowed=1
        intermediate_shell_codes(2)=ket_shell_codes(2)
        intermediate_shell_codes(3)=ket_shell_codes(3)
        intermediate_shell_codes(4)=ket_shell_codes(4)+occupation_change
        intermediate_shell_values(3)=ket_shell_values(3)+half*dble(occupation_change)
        return
    end function select_ket_parity_intermediate_terms
end module jj_intermediate_terms_ket_parity
