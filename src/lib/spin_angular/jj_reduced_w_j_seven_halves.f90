! Look up low-rank pair-tensor elements for j = 7/2 and dispatch higher ranks.
module jj_reduced_w_j_seven_halves
    implicit none (type, external)

    ! Module interface
    private
    public :: rmew7jj
contains
    subroutine rmew7jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero
        use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
            first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
        use jj_reduced_w_j_seven_halves_higher_ranks, only: rmew7bjj
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)  :: bra_term_index, ket_term_index, quasispin_rank, angular_rank
        real(real64), intent(out) :: reduced_w_element

        ! Local declarations
        integer :: w_numerators_0_0(14), w_numerators_1_0(14), w_numerators_0_1(14)

        ! Initial values
        data w_numerators_0_0/32,48,128,80,96,128,20,60,20,108,36,44,156,68/
        data w_numerators_1_0/6,9,120,15,18,24,2*30,0,54,2*0,78,0/
        data w_numerators_0_1/10,35,168,165,286,680,0,30,10,180,60,110,546,408/
        reduced_w_element=zero
        if(first_same_parity_terms_low_j(bra_term_index) /= first_same_parity_terms_low_j(ket_term_index)) return
        if(bra_term_index < 12 .or. bra_term_index > 25) return
        if(quasispin_rank == 0 .and. angular_rank == 0) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_0_0(bra_term_index-11)))
        else if(quasispin_rank == 1 .and. angular_rank == 0) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_1_0(bra_term_index-11)))
        else if(quasispin_rank == 0 .and. angular_rank == 1) then
            if(bra_term_index /= ket_term_index) return
            reduced_w_element=-dsqrt(dble(w_numerators_0_1(bra_term_index-11))/dble(7))
        else
            call rmew7bjj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
        end if
        return
    end subroutine rmew7jj
end module jj_reduced_w_j_seven_halves
