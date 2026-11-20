! Multiply the coupling-chain factors for four ordered active subshells.
module four_subshell_recoupling
    implicit none (type, external)
    private
    public :: reco4
contains
      subroutine reco4(first_shell_position,second_shell_position,third_shell_position,fourth_shell_position, &
          two_rank_1,two_rank_2,two_rank_3,two_second_pair_rank,two_first_pair_rank,evaluate_coefficient, &
          coupling_allowed,recoupling_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: one
      use csf_pair_angular_state,             only: peel_orbital_indices => jlist, bra_subshell_quantum_codes => jjq1, &
          ket_subshell_quantum_codes => jjq2
      use first_active_shell_recoupling, only: diaga1
      use later_active_shell_recoupling, only: diaga2
      use spectator_shell_chain_recoupling, only: diaga3
      use pair_tensor_recoupling, only: diaga4
      implicit none (type, external)
      integer, intent(in)  :: first_shell_position,second_shell_position,third_shell_position, &
          fourth_shell_position,two_rank_1,two_rank_2,two_rank_3,two_second_pair_rank,two_first_pair_rank, &
          evaluate_coefficient
      integer, intent(out) :: coupling_allowed
      real(real64), intent(out) :: recoupling_coefficient
      integer      :: two_bra_shell_j_4, two_ket_shell_j_4, orbital_index_1, orbital_index_2, orbital_index_3, &
          orbital_index_4, coupling_step_count
      real(real64) :: bra_j_weight_product, bra_j_weight_1, bra_j_weight_2, bra_j_weight_3, bra_j_weight_4, &
          diagram_factor
      orbital_index_1=peel_orbital_indices(first_shell_position)
      orbital_index_2=peel_orbital_indices(second_shell_position)
      orbital_index_3=peel_orbital_indices(third_shell_position)
      orbital_index_4=peel_orbital_indices(fourth_shell_position)
      bra_j_weight_1=bra_subshell_quantum_codes(3,orbital_index_1)
      bra_j_weight_2=bra_subshell_quantum_codes(3,orbital_index_2)
      bra_j_weight_3=bra_subshell_quantum_codes(3,orbital_index_3)
      bra_j_weight_4=bra_subshell_quantum_codes(3,orbital_index_4)
      bra_j_weight_product=bra_j_weight_1*bra_j_weight_2*bra_j_weight_3*bra_j_weight_4
      recoupling_coefficient=one/dsqrt(bra_j_weight_product)
      two_bra_shell_j_4=bra_subshell_quantum_codes(3,orbital_index_4)-1
      two_ket_shell_j_4=ket_subshell_quantum_codes(3,orbital_index_4)-1
      recoupling_coefficient=recoupling_coefficient*dsqrt(dble(two_bra_shell_j_4+1))/ &
          dsqrt(dble((two_second_pair_rank+1)*(two_ket_shell_j_4+1)))
      coupling_step_count=third_shell_position-second_shell_position
      if(coupling_step_count > 1) then
        coupling_allowed=0
        call diaga3(second_shell_position,third_shell_position,two_first_pair_rank,evaluate_coefficient, &
            coupling_allowed,diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
      end if
      coupling_step_count=fourth_shell_position-third_shell_position
      if(coupling_step_count > 1) then
        coupling_allowed=0
        call diaga3(third_shell_position,fourth_shell_position,two_second_pair_rank,evaluate_coefficient, &
            coupling_allowed,diagram_factor)
        if(coupling_allowed == 0)return
        recoupling_coefficient=diagram_factor*recoupling_coefficient
      end if
      coupling_allowed=0
      call diaga2(first_shell_position,fourth_shell_position,two_second_pair_rank,evaluate_coefficient, &
          coupling_allowed,diagram_factor)
      if(coupling_allowed == 0)return
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      coupling_allowed=0
      call diaga4(first_shell_position,second_shell_position,two_rank_1,two_rank_2,two_first_pair_rank, &
          evaluate_coefficient,coupling_allowed,diagram_factor)
      if(coupling_allowed == 0)return
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      coupling_allowed=0
      call diaga4(second_shell_position,third_shell_position,two_first_pair_rank,two_rank_3,two_second_pair_rank, &
          evaluate_coefficient,coupling_allowed,diagram_factor)
      if(coupling_allowed == 0)return
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      if(first_shell_position == 1.and.second_shell_position == 2)return
      coupling_allowed=0
      call diaga1(first_shell_position,two_rank_1,evaluate_coefficient,coupling_allowed,diagram_factor)
      if(coupling_allowed == 0)return
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      coupling_step_count=second_shell_position-first_shell_position
      if(first_shell_position == 1)coupling_step_count=second_shell_position-1-first_shell_position
      if(coupling_step_count <= 1)return
      coupling_allowed=0
      call diaga3(first_shell_position,second_shell_position,two_rank_1,evaluate_coefficient,coupling_allowed, &
          diagram_factor)
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      return
      end subroutine reco4
end module four_subshell_recoupling
