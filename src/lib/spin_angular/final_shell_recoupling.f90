! Recouple a tensor through the last subshell of the coupling chain.
module final_shell_recoupling
    implicit none (type, external)
    private
    public :: diaga5
contains
      subroutine diaga5(coupled_shell_count,first_shell_position,two_tensor_rank,evaluate_coefficient, &
          coupling_allowed,recoupling_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero
      use csf_pair_angular_state,             only: peel_orbital_indices => jlist, bra_intermediate_couplings => jjc1, &
          ket_intermediate_couplings => jjc2, bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2
      use wigner_6j_triangle_selection, only: ixjtik
      use wigner_6j, only: sixj
      implicit none (type, external)
      integer, intent(in)       :: coupled_shell_count,first_shell_position,two_tensor_rank,evaluate_coefficient
      integer, intent(inout)      :: coupling_allowed
      real(real64), intent(out) :: recoupling_coefficient
      integer      :: two_bra_coupled_j,two_ket_coupled_j,last_orbital,two_bra_input_j,two_ket_input_j,two_spectator_j
      real(real64) :: diagram_factor
      recoupling_coefficient = zero
      two_bra_coupled_j=bra_intermediate_couplings(coupled_shell_count-1)-1
      two_ket_coupled_j=ket_intermediate_couplings(coupled_shell_count-1)-1
      last_orbital=peel_orbital_indices(coupled_shell_count)
      if(first_shell_position == coupled_shell_count) then
        two_bra_input_j=bra_subshell_quantum_codes(3,last_orbital)-1
        two_ket_input_j=ket_subshell_quantum_codes(3,last_orbital)-1
        two_spectator_j=bra_intermediate_couplings(coupled_shell_count-2)-1
      else
        two_spectator_j=bra_subshell_quantum_codes(3,last_orbital)-1
        two_bra_input_j=bra_intermediate_couplings(coupled_shell_count-2)-1
        two_ket_input_j=ket_intermediate_couplings(coupled_shell_count-2)-1
      end if
      if(evaluate_coefficient == 0) then
        if(ixjtik(two_tensor_rank,two_ket_input_j,two_bra_input_j,two_spectator_j,two_bra_coupled_j, &
            two_ket_coupled_j) /= 0) coupling_allowed=1
      else
        call sixj(two_tensor_rank,two_ket_input_j,two_bra_input_j,two_spectator_j,two_bra_coupled_j, &
            two_ket_coupled_j,0,diagram_factor)
        recoupling_coefficient=diagram_factor*dsqrt(dble((two_bra_input_j+1)*(two_ket_coupled_j+1)))
        if(mod(two_tensor_rank+two_spectator_j+two_ket_input_j+two_bra_coupled_j, &
            4) /= 0)recoupling_coefficient=-recoupling_coefficient
        coupling_allowed=1
        if(first_shell_position == coupled_shell_count)return
        if(mod(two_bra_input_j+two_ket_input_j-two_ket_coupled_j-two_bra_coupled_j+2*two_spectator_j, &
            4) /= 0)recoupling_coefficient=-recoupling_coefficient
      end if
      return
      end subroutine diaga5
end module final_shell_recoupling
