! Form a one-electron tensor coefficient for unchanged occupations at one subshell.
module one_electron_tensor_same_subshell
    implicit none (type, external)
    private
    public :: oneparticlejj1
contains
      subroutine oneparticlejj1(active_shell_count,tensor_rank,bra_csf_index,ket_csf_index,bra_shell_position, &
          ket_shell_position,angular_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, angular_tolerance => eps
      use jj_tensor_coupling_tables, only: ket_shell_values_1 => bd1, ket_shell_values_2 => bd2, ket_shell_values_3 => bd3, &
          ket_shell_values_4 => bd4, bra_shell_values_1 => bk1, bra_shell_values_2 => bk2, &
          bra_shell_values_3 => bk3, bra_shell_values_4 => bk4, ket_shell_codes_1 => id1, ket_shell_codes_2 => id2, &
          ket_shell_codes_3 => id3, ket_shell_codes_4 => id4, bra_shell_codes_1 => ik1, bra_shell_codes_2 => ik2, &
          bra_shell_codes_3 => ik3, bra_shell_codes_4 => ik4
      use one_electron_tensor_recoupling_selection, only: recop00
      use one_electron_same_subshell_recoupling, only: recop1
      use one_electron_two_subshell_recoupling, only: recop2
      use jj_w_matrix_element, only: wj1
      use jj_active_subshell_states, only: perko2
      implicit none (type, external)
      integer, intent(in) :: active_shell_count,tensor_rank,bra_csf_index,ket_csf_index,bra_shell_position, &
          ket_shell_position
      real(real64), intent(out) :: angular_coefficient
      integer      :: coupling_allowed
      real(real64) :: recoupling_factor,w_matrix_element,q_projection_1,q_projection_2
      angular_coefficient=zero
      if(bra_shell_position == ket_shell_position) then
        if(bra_csf_index /= ket_csf_index) then
          call recop00(active_shell_count,bra_shell_position,bra_shell_position,tensor_rank,coupling_allowed)
          if(coupling_allowed == 0)return
        end if
        call recop1(active_shell_count,bra_shell_position,tensor_rank,0,coupling_allowed,recoupling_factor)
        if(coupling_allowed == 0)return
        call perko2(bra_shell_position,bra_shell_position,bra_shell_position,bra_shell_position,1)
        q_projection_1=half
        q_projection_2=-half
        call wj1(bra_shell_codes_1,bra_shell_values_1,ket_shell_codes_1,ket_shell_values_1,tensor_rank, &
            q_projection_1,q_projection_2,w_matrix_element)
        if(dabs(w_matrix_element) > angular_tolerance) then
           call recop1(active_shell_count,bra_shell_position,tensor_rank,1,coupling_allowed,recoupling_factor)
           angular_coefficient=w_matrix_element*recoupling_factor*dsqrt(dble(ket_shell_codes_1(3)+1))/dsqrt(dble(2* &
               tensor_rank+1))
           angular_coefficient=-angular_coefficient
        end if
      end if
      return
      end subroutine oneparticlejj1
end module one_electron_tensor_same_subshell
