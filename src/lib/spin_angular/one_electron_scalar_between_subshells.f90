! Form a scalar one-electron coefficient for transfer between subshells.
module one_electron_scalar_between_subshells
    implicit none (type, external)
    private
    public :: onescalar2
contains
      subroutine onescalar2(bra_csf_index,ket_csf_index,bra_shell_position,ket_shell_position,angular_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, angular_tolerance => eps
      use csf_pair_angular_state, only: bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2, &
          bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2, peel_orbital_indices => jlist, &
          core_orbital_indices => klist, max_orbitals, core_orbital_count => ncore, peel_orbital_count => npeel, &
          bra_occupations => nq1, ket_occupations => nq2
      use jj_tensor_coupling_tables, only: ket_shell_values_1 => bd1, ket_shell_values_2 => bd2, ket_shell_values_3 => bd3, &
          ket_shell_values_4 => bd4, bra_shell_values_1 => bk1, bra_shell_values_2 => bk2, &
          bra_shell_values_3 => bk3, bra_shell_values_4 => bk4, ket_shell_codes_1 => id1, ket_shell_codes_2 => id2, &
          ket_shell_codes_3 => id3, ket_shell_codes_4 => id4, bra_shell_codes_1 => ik1, bra_shell_codes_2 => ik2, &
          bra_shell_codes_3 => ik3, bra_shell_codes_4 => ik4
      use one_electron_scalar_recoupling_selection, only: recoonescalar
      use jj_active_subshell_states, only: perko2
      use two_subshell_recoupling, only: reco2
      use subshell_product_a_a, only: gg12
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_position,ket_shell_position
      real(real64), intent(out) :: angular_coefficient
      integer      :: coupling_allowed,first_shell_position,last_shell_position,intervening_electron_count, &
          last_phase_shell,shell_position,orbital_index,phase_sum
      real(real64) :: q_projection_1,q_projection_2,recoupling_factor,subshell_matrix_element
      angular_coefficient=zero
      if(bra_shell_position == ket_shell_position) return
      if(bra_shell_position < ket_shell_position) then
        first_shell_position=bra_shell_position
        last_shell_position=ket_shell_position
      else
        first_shell_position=ket_shell_position
        last_shell_position=bra_shell_position
      end if
      call recoonescalar(-1,first_shell_position,last_shell_position,last_shell_position,last_shell_position,1, &
          coupling_allowed)
      if(coupling_allowed == 0)return
      q_projection_1=half
      q_projection_2=-half
      call perko2(bra_shell_position,ket_shell_position,bra_shell_position,bra_shell_position,2)
      if(ket_shell_codes_1(3) /= ket_shell_codes_2(3)) return
      call reco2(first_shell_position,last_shell_position,ket_shell_codes_2(3),0,coupling_allowed,recoupling_factor)
      if(coupling_allowed == 0)return
      call gg12(bra_shell_codes_1,bra_shell_codes_2,bra_shell_values_1,bra_shell_values_2,ket_shell_codes_1, &
          ket_shell_codes_2,ket_shell_values_1,ket_shell_values_2,q_projection_1,q_projection_2, &
          subshell_matrix_element)
      if(dabs(subshell_matrix_element) > angular_tolerance) then
         call reco2(first_shell_position,last_shell_position,ket_shell_codes_2(3),1,coupling_allowed, &
             recoupling_factor)
         angular_coefficient=subshell_matrix_element*recoupling_factor*dsqrt(dble(ket_shell_codes_1(3)+1))

         ! The transfer sign depends on the bra occupations between the ordered subshells.
         intervening_electron_count=0
         last_phase_shell=last_shell_position-1
         do shell_position=first_shell_position,last_phase_shell
           orbital_index=peel_orbital_indices(shell_position)
           intervening_electron_count=bra_occupations(orbital_index)+intervening_electron_count
         end do
         if((intervening_electron_count/2)*2 == intervening_electron_count) angular_coefficient=-angular_coefficient
         angular_coefficient=-angular_coefficient
        phase_sum=bra_shell_codes_1(5)*bra_shell_codes_1(4)+bra_shell_codes_2(5)*bra_shell_codes_2(4)- &
            ket_shell_codes_1(5)*ket_shell_codes_1(4)-ket_shell_codes_2(5)*ket_shell_codes_2(4)
        if((phase_sum/4)*4 /= phase_sum)angular_coefficient=-angular_coefficient
      end if
      return
      end subroutine onescalar2
end module one_electron_scalar_between_subshells
