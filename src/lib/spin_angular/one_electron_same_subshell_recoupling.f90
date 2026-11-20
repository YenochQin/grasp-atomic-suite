! Recouple a one-electron tensor acting within one active subshell.
module one_electron_same_subshell_recoupling
    implicit none (type, external)
    private
    public :: recop1
contains
      subroutine recop1(active_shell_count,first_shell_position,tensor_rank,evaluate_coefficient,coupling_allowed, &
          recoupling_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: one
      use csf_pair_angular_state,             only: peel_orbital_indices => jlist, bra_subshell_quantum_codes => jjq1, &
          peel_orbital_count => npeel
      use first_active_shell_recoupling, only: diaga1
      use spectator_shell_chain_recoupling, only: diaga3
      use final_shell_recoupling, only: diaga5
      implicit none (type, external)
      integer, intent(in)       :: active_shell_count, first_shell_position, tensor_rank, evaluate_coefficient
      integer, intent(out)      :: coupling_allowed
      real(real64), intent(out) :: recoupling_coefficient
      integer      :: orbital_index, coupling_step_count, coupled_shell_count
      real(real64) :: bra_subshell_j_weight, diagram_factor
      coupling_allowed=1
      orbital_index=peel_orbital_indices(first_shell_position)
      bra_subshell_j_weight=dble(bra_subshell_quantum_codes(3,orbital_index))
      recoupling_coefficient=one/dsqrt(bra_subshell_j_weight)
      if(peel_orbital_count == 1 .and. active_shell_count == -1)return
      if(active_shell_count == -1) then
         coupled_shell_count = peel_orbital_count
      else
         coupled_shell_count = active_shell_count
      end if
      coupling_allowed=0
      if(evaluate_coefficient /= 0) then
        if(tensor_rank == 0) then
          coupling_allowed=1
          return
        end if
      end if
      coupling_allowed=1
      if(coupled_shell_count == 1) return
      coupling_allowed=0
      if(coupled_shell_count /= 2) then
        call diaga5(coupled_shell_count,first_shell_position,2*tensor_rank,evaluate_coefficient,coupling_allowed, &
            diagram_factor)
        recoupling_coefficient=diagram_factor*recoupling_coefficient
        if(coupling_allowed == 0) return
        if(first_shell_position == coupled_shell_count) return
        coupling_allowed=0
      end if
      call diaga1(first_shell_position,2*tensor_rank,evaluate_coefficient,coupling_allowed,diagram_factor)
      if(coupling_allowed == 0)return
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      if(coupled_shell_count == 2) return
      coupling_step_count=coupled_shell_count-first_shell_position
      if(first_shell_position == 1)coupling_step_count=coupled_shell_count-1-first_shell_position
      if(coupling_step_count <= 1)return
      coupling_allowed=0
      call diaga3(first_shell_position,coupled_shell_count,2*tensor_rank,evaluate_coefficient,coupling_allowed, &
          diagram_factor)
      recoupling_coefficient=diagram_factor*recoupling_coefficient
      return
      end subroutine recop1
end module one_electron_same_subshell_recoupling
