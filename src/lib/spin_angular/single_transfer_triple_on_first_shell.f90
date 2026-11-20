! Evaluate single-electron transfer with three operators acting on subshell 1.
module single_transfer_triple_on_first_shell
    implicit none (type, external)
    private
    public :: el31
contains
      subroutine el31(bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,bra_integral_shell_1, &
          bra_integral_shell_2,ket_integral_shell_1,ket_integral_shell_2,interaction_type)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, angular_tolerance => eps
      use csf_pair_angular_state,             only: &
          bra_occupations => nq1, peel_orbital_indices => jlist, peel_orbital_count => npeel
      use orbital_csf_basis,           only: orbital_kappa
      use jj_tensor_coupling_tables, only: ket_shell_values_1 => bd1, ket_shell_values_2 => bd2, ket_shell_values_3 => bd3, &
          ket_shell_values_4 => bd4, bra_shell_values_1 => bk1, bra_shell_values_2 => bk2, &
          bra_shell_values_3 => bk3, bra_shell_values_4 => bk4, ket_shell_codes_1 => id1, ket_shell_codes_2 => id2, &
          ket_shell_codes_3 => id3, ket_shell_codes_4 => id4, bra_shell_codes_1 => ik1, bra_shell_codes_2 => ik2, &
          bra_shell_codes_3 => ik3, bra_shell_codes_4 => ik4
      use scalar_interaction_recoupling_selection, only: reco
      use two_subshell_recoupling, only: reco2
      use jj_active_subshell_states, only: perko2
      use shared_tensor_rank_range, only: itrexg
      use wigner_6j_triangle_selection, only: ixjtik
      use coulomb_angular_factor, only: coulom
      use breit_radial_rank_ranges, only: snrc
      use wigner_6j, only: sixj
      use mcp_coefficient_output, only: speak
      use subshell_product_a_aw, only: gg1222
      use breit_angular_coefficient_buffer, only: talk
      use breit_radial_integral_coefficients, only: cxk
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,bra_integral_shell_1, &
          bra_integral_shell_2,ket_integral_shell_1,ket_integral_shell_2,interaction_type
      integer ::  orbital_a,orbital_b,work_index,two_multipole_rank_plus_one,two_coupled_rank_plus_one, &
          coupling_allowed,bra_integral_orbital_1,bra_integral_orbital_2,ket_integral_orbital_1, &
          ket_integral_orbital_2,two_rank_range_size,first_coupling_rank_weight,last_coupling_rank_weight, &
                  direct_integral_type,exchange_integral_type,phase_sum,subshell_order_sign,phase_orbital, &
                      first_ordered_shell,last_ordered_shell,phase_end_position,coupled_rank,abs_kappa_1, &
                  abs_kappa_2,multipole_rank,first_direct_nu,direct_nu_count,first_exchange_nu,exchange_nu_count, &
                      breit_rank_index,intervening_electron_count
      integer :: two_orbital_j(2)
      integer :: integral_orbitals(4),two_kappas(4),two_j_plus_one(4)
      real(real64)          :: coupled_term,angular_sum,coulomb_angular_factor,angular_coefficient,six_j_value, &
          q_projection_1,q_projection_2,q_projection_3,q_projection_4,recoupling_factor
      real(real64) :: breit_factors(12)
      if(peel_orbital_count <= 1)return
      bra_integral_orbital_1=peel_orbital_indices(bra_integral_shell_1)
      bra_integral_orbital_2=peel_orbital_indices(bra_integral_shell_2)
      ket_integral_orbital_1=peel_orbital_indices(ket_integral_shell_1)
      ket_integral_orbital_2=peel_orbital_indices(ket_integral_shell_2)
      if(shell_position_a > shell_position_b) then
        first_ordered_shell=shell_position_b
        last_ordered_shell=shell_position_a
      else
        first_ordered_shell=shell_position_a
        last_ordered_shell=shell_position_b
      end if
      call reco(first_ordered_shell,last_ordered_shell,last_ordered_shell,last_ordered_shell,1,coupling_allowed)
      if(coupling_allowed == 0)return
      orbital_a=peel_orbital_indices(shell_position_a)
      orbital_b=peel_orbital_indices(shell_position_b)
      q_projection_1=half
      q_projection_2=half
      q_projection_3=-half
      q_projection_4=-half
      call perko2(shell_position_a,shell_position_b,shell_position_a,shell_position_a,2)
      two_orbital_j(1)=ket_shell_codes_1(3)
      two_orbital_j(2)=ket_shell_codes_2(3)
      abs_kappa_1=(two_orbital_j(1)+1)/2
      abs_kappa_2=(two_orbital_j(2)+1)/2
      call reco2(first_ordered_shell,last_ordered_shell,two_orbital_j(2),0,coupling_allowed,recoupling_factor)
      if(coupling_allowed == 0)return
      first_coupling_rank_weight=itrexg(two_orbital_j(1),two_orbital_j(1),two_orbital_j(1),two_orbital_j(2), &
          two_rank_range_size)+1
      if(two_rank_range_size <= 0)return
      last_coupling_rank_weight=first_coupling_rank_weight+two_rank_range_size-1
      call reco2(first_ordered_shell,last_ordered_shell,two_orbital_j(2),1,coupling_allowed,recoupling_factor)
      if (interaction_type == 2) then
        integral_orbitals(1)=bra_integral_orbital_1
        integral_orbitals(2)=bra_integral_orbital_2
        integral_orbitals(3)=ket_integral_orbital_1
        integral_orbitals(4)=ket_integral_orbital_2
        two_kappas(1)=2*orbital_kappa(integral_orbitals(1))
        two_kappas(2)=2*orbital_kappa(integral_orbitals(2))
        two_kappas(3)=2*orbital_kappa(integral_orbitals(3))
        two_kappas(4)=2*orbital_kappa(integral_orbitals(4))
        two_j_plus_one(1)=iabs(two_kappas(1))
        two_j_plus_one(2)=iabs(two_kappas(2))
        two_j_plus_one(3)=iabs(two_kappas(3))
        two_j_plus_one(4)=iabs(two_kappas(4))
        call snrc(integral_orbitals,two_kappas,two_j_plus_one,first_direct_nu,direct_nu_count,first_exchange_nu, &
            exchange_nu_count,direct_integral_type,exchange_integral_type)
        if(direct_integral_type <= 0)return
      end if
      do two_multipole_rank_plus_one=first_coupling_rank_weight,last_coupling_rank_weight,2
        multipole_rank=(two_multipole_rank_plus_one-1)/2
        if (interaction_type == 1) then
          call coulom(abs_kappa_2,abs_kappa_1,abs_kappa_1,abs_kappa_1,ket_shell_codes_2(5),ket_shell_codes_1(5), &
              ket_shell_codes_1(5),ket_shell_codes_1(5),multipole_rank,coulomb_angular_factor)
          if(dabs(coulomb_angular_factor) < angular_tolerance) cycle
          coulomb_angular_factor=-coulomb_angular_factor
        end if

        ! Sum the 6j-transformed subshell factors over intermediate tensor ranks.
        angular_sum=zero
        do two_coupled_rank_plus_one=first_coupling_rank_weight,last_coupling_rank_weight,2
          coupled_rank=(two_coupled_rank_plus_one-1)/2
          phase_sum=two_orbital_j(2)-coupled_rank+1
          if((phase_sum/2)*2 /= phase_sum) cycle
          if(ixjtik(two_orbital_j(2),two_orbital_j(1),multipole_rank*2,two_orbital_j(1),two_orbital_j(1), &
              coupled_rank*2) == 0)cycle
          call gg1222(bra_shell_codes_2,bra_shell_codes_1,bra_shell_values_2,bra_shell_values_1,ket_shell_codes_2, &
              ket_shell_codes_1,ket_shell_values_2,ket_shell_values_1,coupled_rank,  &
                      q_projection_1,q_projection_2,q_projection_3,q_projection_4,coupled_term)
          if(dabs(coupled_term) < angular_tolerance) cycle
          call sixj(two_orbital_j(2),two_orbital_j(1),multipole_rank*2,two_orbital_j(1),two_orbital_j(1), &
              coupled_rank*2,0,six_j_value)
          coupled_term=coupled_term*six_j_value*dsqrt(dble(two_coupled_rank_plus_one))
          phase_sum=2*two_orbital_j(1)+multipole_rank*2+coupled_rank*2
          if((phase_sum/4)*4 /= phase_sum)coupled_term=-coupled_term
          angular_sum=angular_sum+coupled_term
        end do
        angular_sum=angular_sum*recoupling_factor
        if(dabs(angular_sum) < angular_tolerance) cycle
        subshell_order_sign = 1
        phase_sum=bra_shell_codes_1(5)*bra_shell_codes_1(4)+bra_shell_codes_2(5)*bra_shell_codes_2(4)- &
            ket_shell_codes_1(5)*ket_shell_codes_1(4)-ket_shell_codes_2(5)*ket_shell_codes_2(4)

        ! Apply the orbital-order phase from the occupations and orbital l on each side.
        if((phase_sum/4)*4 /= phase_sum)subshell_order_sign=-subshell_order_sign
        intervening_electron_count=0
        phase_end_position=last_ordered_shell-1
        do work_index=first_ordered_shell,phase_end_position
          phase_orbital=peel_orbital_indices(work_index)
          intervening_electron_count=bra_occupations(phase_orbital)+intervening_electron_count
        end do
        if((intervening_electron_count/2)*2 == intervening_electron_count)angular_sum=-angular_sum
        if (interaction_type == 1) then
           angular_coefficient=coulomb_angular_factor*angular_sum*dble(subshell_order_sign)
           call speak(bra_csf_index,ket_csf_index,bra_integral_orbital_1,bra_integral_orbital_2, &
               ket_integral_orbital_1,ket_integral_orbital_2,multipole_rank,angular_coefficient)
        else if (interaction_type == 2) then
          breit_rank_index=(multipole_rank-first_direct_nu)/2+1
          if(((multipole_rank-first_direct_nu)/2)*2 == (multipole_rank-first_direct_nu)) then
            call cxk(breit_factors,integral_orbitals,two_kappas,multipole_rank,multipole_rank,2,1)
            if(dabs(breit_factors(1)) > angular_tolerance) then
              angular_coefficient=-breit_factors(1)*angular_sum
              if(dabs(angular_coefficient) > angular_tolerance)                                &
              call talk(bra_csf_index,ket_csf_index,multipole_rank,integral_orbitals(1),integral_orbitals(3), &
                  integral_orbitals(2),integral_orbitals(4),3,angular_coefficient)
            end if
          end if
        end if
      end do
      return
      end subroutine el31
end module single_transfer_triple_on_first_shell
