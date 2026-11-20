! Evaluate the 9j factor for two tensors joined at a coupling node.
module pair_tensor_recoupling
    implicit none (type, external)
    private
    public :: diaga4
contains
      subroutine diaga4(first_shell_position,second_shell_position,two_rank_1,two_rank_2,two_tensor_rank, &
          evaluate_coefficient,coupling_allowed,recoupling_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero
      use csf_pair_angular_state,             only: peel_orbital_indices => jlist, bra_subshell_quantum_codes => jjq1, &
          ket_subshell_quantum_codes => jjq2, bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2
      use angular_triangle_and_parity, only: ittk
      use wigner_9j, only: nine
      implicit none (type, external)
      integer, intent(in)       :: first_shell_position,second_shell_position,two_rank_1,two_rank_2, &
          two_tensor_rank,evaluate_coefficient
      integer, intent(inout)      :: coupling_allowed
      real(real64), intent(out) :: recoupling_coefficient
      integer :: orbital_index_1,orbital_index_2,two_bra_shell_j_1,two_bra_shell_j_2,two_ket_shell_j_1, &
          two_ket_shell_j_2,coupled_position,preceding_coupling_position,two_bra_coupled_j,two_ket_coupled_j, &
          two_bra_prefix_j,two_ket_prefix_j
      real(real64) :: nine_j_value
      recoupling_coefficient = zero
      orbital_index_1=peel_orbital_indices(first_shell_position)
      orbital_index_2=peel_orbital_indices(second_shell_position)
      two_bra_shell_j_1=bra_subshell_quantum_codes(3,orbital_index_1)-1
      two_bra_shell_j_2=bra_subshell_quantum_codes(3,orbital_index_2)-1
      two_ket_shell_j_1=ket_subshell_quantum_codes(3,orbital_index_1)-1
      two_ket_shell_j_2=ket_subshell_quantum_codes(3,orbital_index_2)-1
      if(first_shell_position == 1.and.second_shell_position == 2) then
        two_bra_prefix_j=two_bra_shell_j_1
        two_ket_prefix_j=two_ket_shell_j_1
        two_bra_coupled_j=bra_intermediate_couplings(1)-1
        two_ket_coupled_j=ket_intermediate_couplings(1)-1
      else
        coupled_position=second_shell_position-1
        two_bra_coupled_j=bra_intermediate_couplings(coupled_position)-1
        two_ket_coupled_j=ket_intermediate_couplings(coupled_position)-1
        preceding_coupling_position=second_shell_position-2
        two_bra_prefix_j=bra_intermediate_couplings(preceding_coupling_position)-1
        two_ket_prefix_j=ket_intermediate_couplings(preceding_coupling_position)-1
      end if
      if(evaluate_coefficient == 0) then
        coupling_allowed = 0
        if (ittk(two_ket_prefix_j,two_rank_1,two_bra_prefix_j) == 0) return
        if (ittk(two_ket_shell_j_2,two_rank_2,two_bra_shell_j_2) == 0) return
        if (ittk(two_ket_coupled_j,two_tensor_rank,two_bra_coupled_j) == 0) return
        if (ittk(two_ket_prefix_j,two_ket_shell_j_2,two_ket_coupled_j) == 0) return
        if (ittk(two_rank_1,two_rank_2,two_tensor_rank) == 0) return
        if (ittk(two_bra_prefix_j,two_bra_shell_j_2,two_bra_coupled_j) == 0) return
        coupling_allowed = 1
        return
      end if
      call nine(two_ket_prefix_j,two_rank_1,two_bra_prefix_j,two_ket_shell_j_2,two_rank_2,two_bra_shell_j_2, &
          two_ket_coupled_j,two_tensor_rank,two_bra_coupled_j,0,coupling_allowed,nine_j_value)
      recoupling_coefficient=nine_j_value*dsqrt(dble((two_bra_prefix_j+1)*(two_tensor_rank+1)*(two_bra_shell_j_2+1) &
          *(two_ket_coupled_j+1)))
      return
      end subroutine diaga4
end module pair_tensor_recoupling
