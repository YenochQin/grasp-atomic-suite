! Look up pair-tensor elements for j = 3/2 by quasispin and angular rank.
module jj_reduced_w_j_three_halves
    implicit none (type, external)

    ! Module interface
    private
    public :: rmew3jj
contains
    subroutine rmew3jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: error_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
            first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: bra_term_index, ket_term_index, quasispin_rank, angular_rank
        real(real64), intent(out) :: reduced_w_element

        ! Local declarations
        integer               :: bra_table_index
        integer :: w_numerators_0_0(3), w_numerators_1_0(3), w_numerators_0_1(3)

        ! Initial values
        data w_numerators_0_0/16,6,10/
        data w_numerators_1_0/12,12,0/
        data w_numerators_0_1/12,0,12/
        reduced_w_element=zero
        if(first_same_parity_terms_low_j(bra_term_index) /= first_same_parity_terms_low_j(ket_term_index)) return
        if(bra_term_index < 2 .or. bra_term_index > 5) return
        bra_table_index=bra_term_index-2
        if(quasispin_rank == 0 .and. angular_rank == 0) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_0_0(bra_table_index)))
        else if(quasispin_rank == 1 .and. angular_rank == 0) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_1_0(bra_table_index)))
        else if(quasispin_rank == 0 .and. angular_rank == 1) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_0_1(bra_table_index)))
        else if(quasispin_rank == 1 .and. angular_rank == 2) then
            if(bra_term_index == 3 .and. ket_term_index == 3) then
                reduced_w_element=dsqrt(dble(60))
            else if(bra_term_index == 5 .and. ket_term_index == 4) then
                reduced_w_element=dsqrt(dble(30))
            else if(bra_term_index == 4 .and. ket_term_index == 5) then
                reduced_w_element=-dsqrt(dble(30))
            end if
        else if(quasispin_rank == 0 .and. angular_rank == 3) then
            if(bra_term_index == 3 .and. ket_term_index == 3) then
                reduced_w_element=-dsqrt(dble(28))
            else if(bra_term_index == 5 .and. ket_term_index == 5) then
                reduced_w_element=dsqrt(dble(28))
            end if
        else
            write(error_unit,"(A,4I5)") " J1 J2 = ",bra_term_index,ket_term_index
            write(error_unit,"(A)") " ERROR IN SUB. RMEW3JJ "
            stop
        end if
        return
    end subroutine rmew3jj
end module jj_reduced_w_j_three_halves
