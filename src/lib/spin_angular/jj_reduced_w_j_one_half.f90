! Look up pair-tensor elements for j = 1/2 by quasispin and angular rank.
module jj_reduced_w_j_one_half
    implicit none (type, external)
    private
    public :: rmew1jj
contains
      subroutine rmew1jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero
      use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
          first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
      implicit none (type, external)
      integer,      intent(in)  :: bra_term_index, ket_term_index, quasispin_rank, angular_rank
      real(real64), intent(out) :: reduced_w_element
      integer :: w_numerators_1_0(2), w_numerators_0_1(2)
      data w_numerators_0_1/6,0/
      data w_numerators_1_0/0,6/
      reduced_w_element=zero
      if(first_same_parity_terms_low_j(bra_term_index) /= first_same_parity_terms_low_j(ket_term_index)) return
      if(bra_term_index > 2) return
      if(quasispin_rank == 0 .and. angular_rank == 0) then
        reduced_w_element=-dsqrt(dble(2))
      else if(quasispin_rank == 1 .and. angular_rank == 0) then
        reduced_w_element=-dsqrt(dble(w_numerators_1_0(bra_term_index)))
      else if(quasispin_rank == 0 .and. angular_rank == 1) then
        reduced_w_element=-dsqrt(dble(w_numerators_0_1(bra_term_index)))
      else
        write(0,"(A,4I5)") " J1 J2 = ",bra_term_index,ket_term_index
        write(0,"(A)") " ERROR IN SUB. RMEW1JJ "
        stop
      end if
      return
      end subroutine rmew1jj
end module jj_reduced_w_j_one_half
