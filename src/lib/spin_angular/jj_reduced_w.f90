! Select a table of completely reduced pair-tensor matrix elements by orbital j.
module jj_reduced_w
    implicit none (type, external)
    private
    public :: rwjj
contains
      subroutine rwjj(two_orbital_j,bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants, only: eleven, angular_tolerance => eps, four, half, one, seven, ten, tenth, three, two, zero
      use jj_reduced_coefficients, only: last_opposite_parity_terms_low_j => imgnjj, last_same_parity_terms_low_j => imgtjj, &
          first_opposite_parity_terms_low_j => impnjj, first_same_parity_terms_low_j => imptjj
      use jj_reduced_w_j_one_half, only: rmew1jj
      use jj_reduced_w_j_three_halves, only: rmew3jj
      use jj_reduced_w_j_five_halves, only: rmew5jj
      use jj_reduced_w_j_seven_halves, only: rmew7jj
      implicit none (type, external)
      integer,      intent(in)  :: two_orbital_j, bra_term_index, ket_term_index, quasispin_rank, angular_rank
      real(real64), intent(out) :: reduced_w_element
      if(two_orbital_j == 1) then
         call rmew1jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
      else if(two_orbital_j == 3) then
         call rmew3jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
      else if(two_orbital_j == 5) then
         call rmew5jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
      else if(two_orbital_j == 7) then
         call rmew7jj(bra_term_index,ket_term_index,quasispin_rank,angular_rank,reduced_w_element)
      else
         write(0,"(A,I5)") " KLAIDA SUB. RWJJ J=",two_orbital_j
         stop
      end if
      return
      end subroutine rwjj
end module jj_reduced_w
