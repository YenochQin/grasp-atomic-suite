! Evaluate two electrons transferred together between a pair of subshells.
module two_electron_pair_transfer
    implicit none (type, external)
    private
    public :: el2
contains
      subroutine el2(bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,interaction_type)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, half, angular_tolerance => eps
      use csf_pair_angular_state,             only: peel_orbital_indices => jlist, peel_orbital_count => npeel
      use orbital_csf_basis,           only: orbital_kappa
      use jj_tensor_coupling_tables, only: ket_shell_values_1 => bd1, ket_shell_values_2 => bd2, ket_shell_values_3 => bd3, &
          ket_shell_values_4 => bd4, bra_shell_values_1 => bk1, bra_shell_values_2 => bk2, &
          bra_shell_values_3 => bk3, bra_shell_values_4 => bk4, ket_shell_codes_1 => id1, ket_shell_codes_2 => id2, &
          ket_shell_codes_3 => id3, ket_shell_codes_4 => id4, bra_shell_codes_1 => ik1, bra_shell_codes_2 => ik2, &
          bra_shell_codes_3 => ik3, bra_shell_codes_4 => ik4
      use scalar_interaction_recoupling_selection, only: reco
      use two_subshell_recoupling, only: reco2
      use shared_tensor_rank_range, only: itrexg
      use angular_momentum_triangle, only: itrig
      use wigner_6j_triangle_selection, only: ixjtik
      use breit_radial_rank_ranges, only: snrc
      use coulomb_angular_factor, only: coulom
      use subshell_product_w_w, only: gg1122
      use wigner_6j, only: sixj
      use breit_radial_integral_coefficients, only: cxk
      use breit_angular_coefficient_buffer, only: talk
      use jj_active_subshell_states, only: perko2
      use mcp_coefficient_output, only: speak
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,interaction_type
      integer :: coupling_allowed,orbital_a,orbital_b,direct_integral_type,exchange_integral_type, &
          first_coupling_rank_weight,first_interaction_rank_weight,last_coupling_rank_weight, &
          last_interaction_rank_weight,two_rank_range_size,work_index,two_multipole_rank_plus_one, &
          two_coupled_rank_plus_one, &
                 phase_sum,orbital_phase_sum,subshell_order_sign,coupled_rank,first_ordered_shell, &
                     last_ordered_shell, &
                 multipole_rank,abs_kappa_1,abs_kappa_2,breit_rank_index,first_direct_nu,direct_nu_count, &
                     first_exchange_nu,exchange_nu_count,next_breit_radial_rank,breit_radial_rank,breit_component
      integer :: two_orbital_j(2)
      integer :: integral_orbitals(4),two_kappas(4),two_j_plus_one(4)
      real(real64)          :: coupled_term,coulomb_angular_factor,angular_sum,angular_coefficient,q_projection_1, &
          q_projection_2,q_projection_3,q_projection_4,recoupling_factor,six_j_value
      real(real64)    :: breit_factors(12)
      real(real64) :: direct_breit_coefficients(12, 20)
      if(peel_orbital_count <= 1)return
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
      first_coupling_rank_weight=itrexg(two_orbital_j(1),two_orbital_j(1),two_orbital_j(2),two_orbital_j(2), &
          two_rank_range_size)+1
      if(two_rank_range_size <= 0)return
      last_coupling_rank_weight=first_coupling_rank_weight+two_rank_range_size-1
      first_interaction_rank_weight=itrexg(two_orbital_j(1),two_orbital_j(2),two_orbital_j(1),two_orbital_j(2), &
          two_rank_range_size)+1
      if(two_rank_range_size <= 0) return
      last_interaction_rank_weight=first_interaction_rank_weight+two_rank_range_size-1
      if (interaction_type == 2) then
        integral_orbitals(1)=orbital_a
        integral_orbitals(2)=orbital_a
        integral_orbitals(3)=orbital_b
        integral_orbitals(4)=orbital_b
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
        do work_index=1,20
          direct_breit_coefficients(1,work_index) =zero
          direct_breit_coefficients(2,work_index) =zero
          direct_breit_coefficients(3,work_index) =zero
          direct_breit_coefficients(4,work_index) =zero
          direct_breit_coefficients(5,work_index) =zero
          direct_breit_coefficients(6,work_index) =zero
          direct_breit_coefficients(7,work_index) =zero
          direct_breit_coefficients(8,work_index) =zero
          direct_breit_coefficients(9,work_index) =zero
          direct_breit_coefficients(10,work_index)=zero
          direct_breit_coefficients(11,work_index)=zero
          direct_breit_coefficients(12,work_index)=zero
        end do
      end if
      do two_multipole_rank_plus_one=first_interaction_rank_weight,last_interaction_rank_weight,2
        multipole_rank=(two_multipole_rank_plus_one-1)/2
        if (interaction_type == 1) then
          call coulom(abs_kappa_1,abs_kappa_1,abs_kappa_2,abs_kappa_2,ket_shell_codes_1(5),ket_shell_codes_1(5), &
              ket_shell_codes_2(5),ket_shell_codes_2(5),multipole_rank,coulomb_angular_factor)
          if(dabs(coulomb_angular_factor) < angular_tolerance) cycle
          coulomb_angular_factor=-half*coulomb_angular_factor
        end if

        ! Sum the 6j-transformed subshell factors over intermediate tensor ranks.
        angular_sum=zero
          do two_coupled_rank_plus_one=first_coupling_rank_weight,last_coupling_rank_weight,2
            coupled_rank=(two_coupled_rank_plus_one-1)/2
            call reco2(first_ordered_shell,last_ordered_shell,coupled_rank*2,0,coupling_allowed,recoupling_factor)
            if(coupling_allowed /= 0) then
              if(ixjtik(two_orbital_j(1),two_orbital_j(2),multipole_rank*2,two_orbital_j(2),two_orbital_j(1), &
                  coupled_rank*2) /= 0) then
                call gg1122(coupled_rank,coupled_rank,q_projection_1,q_projection_2,q_projection_3,q_projection_4, &
                    coupled_term)
                if(dabs(coupled_term) > angular_tolerance) then
                  call reco2(first_ordered_shell,last_ordered_shell,coupled_rank*2,1,coupling_allowed, &
                      recoupling_factor)
                  coupled_term=coupled_term*recoupling_factor
                  call sixj(two_orbital_j(1),two_orbital_j(2),multipole_rank*2,two_orbital_j(2),two_orbital_j(1), &
                      coupled_rank*2,0,six_j_value)
                  coupled_term=coupled_term*six_j_value*dsqrt(dble(two_coupled_rank_plus_one))
                  phase_sum=bra_shell_codes_1(3)+bra_shell_codes_2(3)+multipole_rank*2+coupled_rank*2
                  if((phase_sum/4)*4 /= phase_sum)coupled_term=-coupled_term
                  angular_sum=angular_sum+coupled_term
                end if
              end if
            end if
          end do
        subshell_order_sign=1
        orbital_phase_sum=bra_shell_codes_1(5)*bra_shell_codes_1(4)+bra_shell_codes_2(5)*bra_shell_codes_2(4)- &
            ket_shell_codes_1(5)*ket_shell_codes_1(4)-ket_shell_codes_2(5)*ket_shell_codes_2(4)

        ! Apply the orbital-order phase from the occupations and orbital l on each side.
        if((orbital_phase_sum/4)*4 /= orbital_phase_sum)subshell_order_sign=-subshell_order_sign
        if (interaction_type == 1) then
          angular_coefficient=coulomb_angular_factor*angular_sum*dble(subshell_order_sign)
          if(dabs(angular_coefficient) > angular_tolerance)call speak(bra_csf_index,ket_csf_index,orbital_a, &
              orbital_a,orbital_b,orbital_b,multipole_rank,angular_coefficient)
        else if (interaction_type == 2) then
          breit_radial_rank=multipole_rank
          if(((breit_radial_rank-first_direct_nu)/2)*2 == (breit_radial_rank-first_direct_nu)) then
            if((itrig(two_j_plus_one(1),two_j_plus_one(3),breit_radial_rank+breit_radial_rank+1) /= 0) .and.   &
               (itrig(two_j_plus_one(2),two_j_plus_one(4),breit_radial_rank+breit_radial_rank+1) /= 0)) then
              breit_rank_index=(breit_radial_rank-first_direct_nu)/2+1
              if(breit_radial_rank > 0) then
                call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,1)
                do breit_component = 1,4
                  direct_breit_coefficients(breit_component,breit_rank_index)= &
                      direct_breit_coefficients(breit_component,breit_rank_index)-half*angular_sum* &
                      breit_factors(breit_component)
                end do
              end if
            end if
          end if
          breit_radial_rank=multipole_rank+1
          if(((breit_radial_rank-first_direct_nu)/2)*2 == (breit_radial_rank-first_direct_nu)) then
            if((itrig(two_j_plus_one(1),two_j_plus_one(3),breit_radial_rank+breit_radial_rank-1) /= 0) .and.  &
               (itrig(two_j_plus_one(2),two_j_plus_one(4),breit_radial_rank+breit_radial_rank-1) /= 0)) then
              breit_rank_index=(breit_radial_rank-first_direct_nu)/2+1
              if(breit_rank_index <= direct_nu_count) then
                call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,1)
                do breit_component = 1,4
                  direct_breit_coefficients(breit_component,breit_rank_index)= &
                      direct_breit_coefficients(breit_component,breit_rank_index)-half*angular_sum* &
                      breit_factors(breit_component)
                end do
              end if
            end if
          end if
          breit_radial_rank=multipole_rank-1
          if(((breit_radial_rank-first_direct_nu)/2)*2 == (breit_radial_rank-first_direct_nu)) then
            if((itrig(two_j_plus_one(1),two_j_plus_one(3),breit_radial_rank+breit_radial_rank+3) /= 0) .and.   &
               (itrig(two_j_plus_one(2),two_j_plus_one(4),breit_radial_rank+breit_radial_rank+3) /= 0)) then
              if(breit_radial_rank >= 0) then
                breit_rank_index=(breit_radial_rank-first_direct_nu)/2+1
                if(breit_rank_index < direct_nu_count) then
                  call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,1)
                  do breit_component = 1,12
                    direct_breit_coefficients(breit_component,breit_rank_index)= &
                        direct_breit_coefficients(breit_component,breit_rank_index)-half*angular_sum* &
                        breit_factors(breit_component)
                  end do
                end if
              end if
            end if
          end if
        end if
      end do
      if (interaction_type == 2) then
        do breit_rank_index = 1,direct_nu_count
          breit_radial_rank=first_direct_nu+2*(breit_rank_index-1)
          call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(1),integral_orbitals(3), &
              integral_orbitals(2),integral_orbitals(4),1,direct_breit_coefficients(1,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(3),integral_orbitals(1), &
              integral_orbitals(4),integral_orbitals(2),1,direct_breit_coefficients(2,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(1),integral_orbitals(3), &
              integral_orbitals(4),integral_orbitals(2),1,direct_breit_coefficients(3,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(3),integral_orbitals(1), &
              integral_orbitals(2),integral_orbitals(4),1,direct_breit_coefficients(4,breit_rank_index))
          if(breit_rank_index == direct_nu_count) cycle
          next_breit_radial_rank=breit_radial_rank+1
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(1),integral_orbitals(3), &
              integral_orbitals(2),integral_orbitals(4),2,direct_breit_coefficients(5,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(2),integral_orbitals(4), &
              integral_orbitals(1),integral_orbitals(3),2,direct_breit_coefficients(6,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(3),integral_orbitals(1), &
              integral_orbitals(4),integral_orbitals(2),2,direct_breit_coefficients(7,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(4),integral_orbitals(2), &
              integral_orbitals(3),integral_orbitals(1),2,direct_breit_coefficients(8,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(1),integral_orbitals(3), &
              integral_orbitals(4),integral_orbitals(2),2,direct_breit_coefficients(9,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(4),integral_orbitals(2), &
              integral_orbitals(1),integral_orbitals(3),2,direct_breit_coefficients(10,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(3),integral_orbitals(1), &
              integral_orbitals(2),integral_orbitals(4),2,direct_breit_coefficients(11,breit_rank_index))
          call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(2),integral_orbitals(4), &
              integral_orbitals(3),integral_orbitals(1),2,direct_breit_coefficients(12,breit_rank_index))
        end do
      end if
      return
      end subroutine el2
end module two_electron_pair_transfer
