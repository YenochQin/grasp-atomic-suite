! Reduce the one-to-two-electron tensor matrix element in quasispin.
module jj_reduced_cfp_two_electron_terms
    implicit none (type, external)
    private
    public :: rmeajj11
contains
      subroutine rmeajj11(two_electron_term_index,one_electron_term_index,two_orbital_j,reduced_cfp)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, angular_tolerance => eps
      use jj_subshell_term_decoding, only: rumtjj
      use quasispin_cg_rank_half, only: c0t5s
      implicit none (type, external)
      integer,      intent(in)  :: two_orbital_j, two_electron_term_index, one_electron_term_index
      real(real64), intent(out) :: reduced_cfp
      integer      :: two_bra_q, two_bra_seniority, two_ket_q, two_ket_seniority, two_bra_term_j, two_ket_term_j, &
          bra_occupation, phase_exponent
      real(real64) :: normalization_squared, cg_coefficient, bra_q, operator_q_projection, ket_q, bra_q_projection, &
          ket_q_projection
      reduced_cfp=zero
      call rumtjj(two_electron_term_index,two_orbital_j,two_bra_q,two_bra_seniority,two_bra_term_j)
      call rumtjj(one_electron_term_index,two_orbital_j,two_ket_q,two_ket_seniority,two_ket_term_j)
      operator_q_projection=half
      bra_occupation=2
      bra_q=half*dble(two_bra_q)
      ket_q=half*dble(two_ket_q)
      bra_q_projection=-half*dble(half*(two_orbital_j+1)-bra_occupation)
      ket_q_projection=-half*dble(half*(two_orbital_j+1)-(bra_occupation-1))
      call c0t5s(ket_q,ket_q_projection,operator_q_projection,bra_q,bra_q_projection,cg_coefficient)
      if(dabs(cg_coefficient) < angular_tolerance) return
      normalization_squared=dble(bra_occupation*(two_bra_q+1)*(two_bra_term_j+1))
      reduced_cfp=dsqrt(normalization_squared)/cg_coefficient
      phase_exponent=-bra_occupation-1
      if((phase_exponent/2)*2 /= phase_exponent)reduced_cfp=-reduced_cfp
      return
      end subroutine rmeajj11
end module jj_reduced_cfp_two_electron_terms
