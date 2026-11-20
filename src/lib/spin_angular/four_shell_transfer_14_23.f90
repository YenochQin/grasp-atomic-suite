! Evaluate the four-subshell two-electron transfer with pairs (1,4) and (2,3).
module four_shell_transfer_14_23
    implicit none (type, external)

    ! Module interface
    private
    public :: el53
contains
    subroutine el53(bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,shell_position_c, &
        shell_position_d,orbital_order_case,interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
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
        use four_subshell_recoupling, only: reco4
        use jj_active_subshell_states, only: perko2
        use breit_radial_rank_ranges, only: snrc
        use subshell_product_a_a_a_a, only: gg1234
        use shared_tensor_rank_range, only: itrexg
        use wigner_6j_triangle_selection, only: ixjtik
        use coulomb_angular_factor, only: coulom
        use wigner_6j, only: sixj
        use mcp_coefficient_output, only: speak
        use angular_momentum_triangle, only: itrig
        use breit_radial_integral_coefficients, only: cxk
        use breit_angular_coefficient_buffer, only: talk
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,shell_position_c, &
            shell_position_d,orbital_order_case,interaction_type

        ! Local declarations
        integer :: orbital_a,orbital_b,orbital_c,orbital_d,direct_integral_type,exchange_integral_type,work_index, &
            first_coupling_rank_weight,first_interaction_rank_weight,last_coupling_rank_weight, &
            last_interaction_rank_weight,two_rank_range_size,two_multipole_rank_plus_one,two_coupled_rank_plus_one,  &
            two_table_rank_plus_one,phase_sum,permutation_sign,subshell_order_sign,phase_orbital, &
            coupling_allowed,multipole_rank,rank_table_position,abs_kappa_1,abs_kappa_2,abs_kappa_3, &
            abs_kappa_4, &
            coupled_rank,phase_end_before_b,phase_end_before_d,first_direct_nu,direct_nu_count, &
            first_exchange_nu,exchange_nu_count,breit_rank_index,intervening_electron_count, &
            breit_radial_rank,next_breit_radial_rank,breit_component
        integer :: two_orbital_j(4)
        integer :: integral_orbitals(4),two_kappas(4),two_j_plus_one(4)
        real(real64)          :: coupled_term,angular_sum,coulomb_angular_factor,angular_coefficient,q_projection_1, &
            q_projection_2,q_projection_3,q_projection_4,subshell_matrix_element,recoupling_factor,six_j_value
        real(real64) :: breit_factors(12)
        real(real64) :: recoupling_by_rank(30)
        real(real64) :: direct_breit_coefficients(12, 20),exchange_breit_coefficients(12, 20)
        if(peel_orbital_count <= 3)return
        call reco(shell_position_a,shell_position_d,shell_position_c,shell_position_b,3,coupling_allowed)
        if(coupling_allowed == 0)return
        orbital_a=peel_orbital_indices(shell_position_a)
        orbital_b=peel_orbital_indices(shell_position_b)
        orbital_c=peel_orbital_indices(shell_position_c)
        orbital_d=peel_orbital_indices(shell_position_d)
        if(orbital_order_case == 2) then
            q_projection_1=-half
            q_projection_2=half
            q_projection_3=half
            q_projection_4=-half
        else
            q_projection_1=half
            q_projection_2=-half
            q_projection_3=-half
            q_projection_4=half
        end if
        call perko2(shell_position_a,shell_position_b,shell_position_c,shell_position_d,4)
        two_orbital_j(1)=ket_shell_codes_1(3)
        two_orbital_j(2)=ket_shell_codes_2(3)
        two_orbital_j(3)=ket_shell_codes_3(3)
        two_orbital_j(4)=ket_shell_codes_4(3)
        abs_kappa_1=(two_orbital_j(1)+1)/2
        abs_kappa_2=(two_orbital_j(2)+1)/2
        abs_kappa_3=(two_orbital_j(3)+1)/2
        abs_kappa_4=(two_orbital_j(4)+1)/2
        if (interaction_type == 2) then
            if(orbital_order_case == 1) then
                integral_orbitals(1)=orbital_a
                integral_orbitals(2)=orbital_d
                integral_orbitals(3)=orbital_b
                integral_orbitals(4)=orbital_c
            else if(orbital_order_case == 2) then
                integral_orbitals(1)=orbital_b
                integral_orbitals(2)=orbital_c
                integral_orbitals(3)=orbital_a
                integral_orbitals(4)=orbital_d
            end if
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
            if(direct_integral_type <= 0 .and. exchange_integral_type <= 0)return
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
                exchange_breit_coefficients(1,work_index) =zero
                exchange_breit_coefficients(2,work_index) =zero
                exchange_breit_coefficients(3,work_index) =zero
                exchange_breit_coefficients(4,work_index) =zero
                exchange_breit_coefficients(5,work_index) =zero
                exchange_breit_coefficients(6,work_index) =zero
                exchange_breit_coefficients(7,work_index) =zero
                exchange_breit_coefficients(8,work_index) =zero
                exchange_breit_coefficients(9,work_index) =zero
                exchange_breit_coefficients(10,work_index)=zero
                exchange_breit_coefficients(11,work_index)=zero
                exchange_breit_coefficients(12,work_index)=zero
            end do
        end if
        call gg1234(bra_shell_codes_1,bra_shell_codes_2,bra_shell_codes_3,bra_shell_codes_4,bra_shell_values_1, &
            bra_shell_values_2,bra_shell_values_3,bra_shell_values_4,ket_shell_codes_1,ket_shell_codes_2,  &
            ket_shell_codes_3,ket_shell_codes_4,ket_shell_values_1,ket_shell_values_2,ket_shell_values_3, &
            ket_shell_values_4,q_projection_1,q_projection_2,q_projection_3,q_projection_4,subshell_matrix_element)
        if(dabs(subshell_matrix_element) < angular_tolerance) return
        first_coupling_rank_weight=itrexg(two_orbital_j(1),two_orbital_j(2),two_orbital_j(3),two_orbital_j(4), &
            two_rank_range_size)+1
        if(two_rank_range_size <= 0)return
        last_coupling_rank_weight=first_coupling_rank_weight+two_rank_range_size-1
        do two_table_rank_plus_one=first_coupling_rank_weight,last_coupling_rank_weight,2
            multipole_rank=(two_table_rank_plus_one-1)/2
            rank_table_position=multipole_rank+1
            if(rank_table_position > 30) call report_rank_limit()

            ! Cache the subshell recoupling factor at each intermediate rank for later orbital pairings.
            recoupling_by_rank(rank_table_position)=zero
            call reco4(shell_position_a,shell_position_b,shell_position_c,shell_position_d,two_orbital_j(1), &
                two_orbital_j(2),two_orbital_j(3),two_orbital_j(4),multipole_rank*2,0,coupling_allowed,recoupling_factor)
            if(coupling_allowed == 0) cycle
            call reco4(shell_position_a,shell_position_b,shell_position_c,shell_position_d,two_orbital_j(1), &
                two_orbital_j(2),two_orbital_j(3),two_orbital_j(4),multipole_rank*2,1,coupling_allowed,recoupling_factor)
            recoupling_by_rank(rank_table_position)=recoupling_factor
        end do
        subshell_order_sign=1
        phase_sum=bra_shell_codes_1(5)*bra_shell_codes_1(4)+bra_shell_codes_2(5)*bra_shell_codes_2(4)- &
            ket_shell_codes_1(5)*ket_shell_codes_1(4)-ket_shell_codes_2(5)*ket_shell_codes_2(4) &
            +bra_shell_codes_3(5)*bra_shell_codes_3(4)-ket_shell_codes_3(5)*ket_shell_codes_3(4)+bra_shell_codes_4(5)* &
            bra_shell_codes_4(4)-ket_shell_codes_4(5)*ket_shell_codes_4(4)

        ! Apply the orbital-order phase from the occupations and orbital l on each side.
        if((phase_sum/4)*4 /= phase_sum)subshell_order_sign=-subshell_order_sign
        intervening_electron_count=0
        phase_end_before_b=shell_position_b-1
        ! The fermionic phase uses the occupation sums over [a,b) and [c,d) in peel order.
        permutation_sign=1
        do work_index=shell_position_a,phase_end_before_b
            phase_orbital=peel_orbital_indices(work_index)
            intervening_electron_count=bra_occupations(phase_orbital)+intervening_electron_count
        end do
        if((intervening_electron_count/2)*2 == intervening_electron_count)permutation_sign=-permutation_sign
        intervening_electron_count=0
        phase_end_before_d=shell_position_d-1
        do work_index=shell_position_c,phase_end_before_d
            phase_orbital=peel_orbital_indices(work_index)
            intervening_electron_count=bra_occupations(phase_orbital)+intervening_electron_count
        end do
        if((intervening_electron_count/2)*2 == intervening_electron_count)permutation_sign=-permutation_sign
        do two_coupled_rank_plus_one=first_coupling_rank_weight,last_coupling_rank_weight,2
            multipole_rank=(two_coupled_rank_plus_one-1)/2
            if (interaction_type == 1) then
                if(orbital_order_case == 2) then
                    call coulom(abs_kappa_2,abs_kappa_3,abs_kappa_1,abs_kappa_4,ket_shell_codes_2(5),ket_shell_codes_3(5), &
                        ket_shell_codes_1(5),ket_shell_codes_4(5),multipole_rank,coulomb_angular_factor)
                else
                    call coulom(abs_kappa_1,abs_kappa_4,abs_kappa_2,abs_kappa_3,ket_shell_codes_1(5),ket_shell_codes_4(5), &
                        ket_shell_codes_2(5),ket_shell_codes_3(5),multipole_rank,coulomb_angular_factor)
                end if
                if(dabs(coulomb_angular_factor) < angular_tolerance) cycle
            end if
            rank_table_position=multipole_rank+1
            if(rank_table_position > 30) call report_rank_limit()
            coupled_term=recoupling_by_rank(rank_table_position)
            if(dabs(coupled_term) < angular_tolerance) cycle
            coupled_term=coupled_term*subshell_matrix_element
            if(dabs(coupled_term) < angular_tolerance) cycle
            coupled_term=coupled_term/dsqrt(dble(two_coupled_rank_plus_one))
            phase_sum=two_orbital_j(4)+two_orbital_j(3)-2*multipole_rank+2
            if(orbital_order_case == 2)phase_sum=two_orbital_j(1)+two_orbital_j(2)-2*multipole_rank+2
            if((phase_sum/4)*4 /= phase_sum)coupled_term=-coupled_term
            angular_sum=coupled_term*dble(permutation_sign)
            if (interaction_type == 1) then
                angular_coefficient=coulomb_angular_factor*angular_sum*dble(subshell_order_sign)
                if(orbital_order_case == 1)call speak(bra_csf_index,ket_csf_index,orbital_a,orbital_d,orbital_b, &
                    orbital_c,multipole_rank,angular_coefficient)
                if(orbital_order_case == 2)call speak(bra_csf_index,ket_csf_index,orbital_b,orbital_c,orbital_a, &
                    orbital_d,multipole_rank,angular_coefficient)
            else if (interaction_type == 2) then
                breit_radial_rank=multipole_rank
                if(((breit_radial_rank-first_direct_nu)/2)*2 == (breit_radial_rank-first_direct_nu)) then
                    if((itrig(two_j_plus_one(1),two_j_plus_one(3),breit_radial_rank+breit_radial_rank+1) /= 0) .and. &
                        (itrig(two_j_plus_one(2),two_j_plus_one(4),breit_radial_rank+breit_radial_rank+1) /= 0)) then
                        breit_rank_index=(breit_radial_rank-first_direct_nu)/2+1
                        if(breit_radial_rank > 0) then
                            call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,1)
                            do breit_component = 1,4
                                direct_breit_coefficients(breit_component,breit_rank_index)= &
                                    direct_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
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
                                    direct_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
                                    breit_factors(breit_component)
                            end do
                        end if
                    end if
                end if
                breit_radial_rank=multipole_rank-1
                if(((breit_radial_rank-first_direct_nu)/2)*2 == (breit_radial_rank-first_direct_nu)) then
                    if((itrig(two_j_plus_one(1),two_j_plus_one(3),breit_radial_rank+breit_radial_rank+3) /= 0) .and.  &
                        (itrig(two_j_plus_one(2),two_j_plus_one(4),breit_radial_rank+breit_radial_rank+3) /= 0)) then
                        if(breit_radial_rank >= 0) then
                            breit_rank_index=(breit_radial_rank-first_direct_nu)/2+1
                            if(breit_rank_index < direct_nu_count) then
                                call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,1)
                                do breit_component = 1,12
                                    direct_breit_coefficients(breit_component,breit_rank_index)= &
                                        direct_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
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
        first_interaction_rank_weight=itrexg(two_orbital_j(1),two_orbital_j(3),two_orbital_j(4),two_orbital_j(2), &
            two_rank_range_size)+1
        if(two_rank_range_size <= 0) return
        last_interaction_rank_weight=first_interaction_rank_weight+two_rank_range_size-1
        do two_multipole_rank_plus_one=first_interaction_rank_weight,last_interaction_rank_weight,2
            multipole_rank=(two_multipole_rank_plus_one-1)/2
            if (interaction_type == 1) then
                if(orbital_order_case == 2) then
                    call coulom(abs_kappa_2,abs_kappa_3,abs_kappa_4,abs_kappa_1,ket_shell_codes_2(5),ket_shell_codes_3(5), &
                        ket_shell_codes_4(5),ket_shell_codes_1(5),multipole_rank,coulomb_angular_factor)
                else
                    call coulom(abs_kappa_1,abs_kappa_4,abs_kappa_3,abs_kappa_2,ket_shell_codes_1(5),ket_shell_codes_4(5), &
                        ket_shell_codes_3(5),ket_shell_codes_2(5),multipole_rank,coulomb_angular_factor)
                end if
                if(dabs(coulomb_angular_factor) < angular_tolerance) cycle
            end if

            ! Sum the 6j-transformed subshell factors over intermediate tensor ranks.
            angular_sum=zero
            do two_coupled_rank_plus_one=first_coupling_rank_weight,last_coupling_rank_weight,2
                coupled_rank=(two_coupled_rank_plus_one-1)/2
                phase_sum=two_orbital_j(1)+two_orbital_j(4)+2*coupled_rank
                if(orbital_order_case == 2)phase_sum=two_orbital_j(2)+two_orbital_j(3)+2*coupled_rank
                if((phase_sum/2)*2 /= phase_sum) cycle
                rank_table_position=coupled_rank+1
                if(rank_table_position > 30) call report_rank_limit()
                coupled_term=recoupling_by_rank(rank_table_position)
                if(dabs(coupled_term) < angular_tolerance) cycle
                coupled_term=coupled_term*subshell_matrix_element
                if(dabs(coupled_term) < angular_tolerance) cycle
                if(ixjtik(two_orbital_j(1),two_orbital_j(3),multipole_rank*2,two_orbital_j(4),two_orbital_j(2), &
                    coupled_rank*2) == 0) cycle
                call sixj(two_orbital_j(1),two_orbital_j(3),multipole_rank*2,two_orbital_j(4),two_orbital_j(2), &
                    coupled_rank*2,0,six_j_value)
                coupled_term=coupled_term*six_j_value*dsqrt(dble(two_coupled_rank_plus_one))
                phase_sum=two_orbital_j(3)-two_orbital_j(4)-4*multipole_rank+2*coupled_rank
                if(orbital_order_case == 2)phase_sum=two_orbital_j(1)-two_orbital_j(2)-4*multipole_rank-2*coupled_rank
                if((phase_sum/4)*4 /= phase_sum)coupled_term=-coupled_term
                angular_sum=angular_sum+coupled_term
            end do
            if(dabs(angular_sum) < angular_tolerance) cycle
            angular_sum=angular_sum*dble(permutation_sign)
            if (interaction_type == 1) then
                angular_coefficient=coulomb_angular_factor*angular_sum*dble(subshell_order_sign)
                if(orbital_order_case == 1)call speak(bra_csf_index,ket_csf_index,orbital_a,orbital_d,orbital_c, &
                    orbital_b,multipole_rank,angular_coefficient)
                if(orbital_order_case == 2)call speak(bra_csf_index,ket_csf_index,orbital_b,orbital_c,orbital_d, &
                    orbital_a,multipole_rank,angular_coefficient)
            else if (interaction_type == 2) then
                breit_radial_rank=multipole_rank
                if(((breit_radial_rank-first_exchange_nu)/2)*2 == (breit_radial_rank-first_exchange_nu)) then
                    if((itrig(two_j_plus_one(1),two_j_plus_one(4),breit_radial_rank+breit_radial_rank+1) /= 0) .and.  &
                        (itrig(two_j_plus_one(2),two_j_plus_one(3),breit_radial_rank+breit_radial_rank+1) /= 0)) then
                        breit_rank_index=(breit_radial_rank-first_exchange_nu)/2+1
                        if(breit_radial_rank > 0) then
                            call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,2)
                            do breit_component = 1,4
                                exchange_breit_coefficients(breit_component,breit_rank_index)= &
                                    exchange_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
                                    breit_factors(breit_component)
                            end do
                        end if
                    end if
                end if
                breit_radial_rank=multipole_rank+1
                if(((breit_radial_rank-first_exchange_nu)/2)*2 == (breit_radial_rank-first_exchange_nu)) then
                    if((itrig(two_j_plus_one(1),two_j_plus_one(4),breit_radial_rank+breit_radial_rank-1) /= 0) .and.  &
                        (itrig(two_j_plus_one(2),two_j_plus_one(3),breit_radial_rank+breit_radial_rank-1) /= 0)) then
                        breit_rank_index=(breit_radial_rank-first_exchange_nu)/2+1
                        if(breit_rank_index <= exchange_nu_count) then
                            call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,2)
                            do breit_component = 1,4
                                exchange_breit_coefficients(breit_component,breit_rank_index)= &
                                    exchange_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
                                    breit_factors(breit_component)
                            end do
                        end if
                    end if
                end if
                breit_radial_rank=multipole_rank-1
                if(((breit_radial_rank-first_exchange_nu)/2)*2 == (breit_radial_rank-first_exchange_nu)) then
                    if((itrig(two_j_plus_one(1),two_j_plus_one(4),breit_radial_rank+breit_radial_rank+3) /= 0) .and.  &
                        (itrig(two_j_plus_one(2),two_j_plus_one(3),breit_radial_rank+breit_radial_rank+3) /= 0)) then
                        if(breit_radial_rank >= 0) then
                            breit_rank_index=(breit_radial_rank-first_exchange_nu)/2+1
                            if(breit_rank_index < exchange_nu_count) then
                                call cxk(breit_factors,integral_orbitals,two_kappas,breit_radial_rank,multipole_rank,1,2)
                                do breit_component = 1,12
                                    exchange_breit_coefficients(breit_component,breit_rank_index)= &
                                        exchange_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
                                        breit_factors(breit_component)
                                end do
                            end if
                        end if
                    end if
                end if
            end if
        end do
        if (interaction_type == 2) then
            do breit_rank_index = 1,exchange_nu_count
                breit_radial_rank=first_exchange_nu+2*(breit_rank_index-1)
                call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(1),integral_orbitals(4), &
                    integral_orbitals(2),integral_orbitals(3),1,exchange_breit_coefficients(1,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(4),integral_orbitals(1), &
                    integral_orbitals(3),integral_orbitals(2),1,exchange_breit_coefficients(2,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(1),integral_orbitals(4), &
                    integral_orbitals(3),integral_orbitals(2),1,exchange_breit_coefficients(3,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,breit_radial_rank,integral_orbitals(4),integral_orbitals(1), &
                    integral_orbitals(2),integral_orbitals(3),1,exchange_breit_coefficients(4,breit_rank_index))
                if(breit_rank_index == exchange_nu_count) cycle
                next_breit_radial_rank=breit_radial_rank+1
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(1),integral_orbitals(4), &
                    integral_orbitals(2),integral_orbitals(3),2,exchange_breit_coefficients(5,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(2),integral_orbitals(3), &
                    integral_orbitals(1),integral_orbitals(4),2,exchange_breit_coefficients(6,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(4),integral_orbitals(1), &
                    integral_orbitals(3),integral_orbitals(2),2,exchange_breit_coefficients(7,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(3),integral_orbitals(2), &
                    integral_orbitals(4),integral_orbitals(1),2,exchange_breit_coefficients(8,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(1),integral_orbitals(4), &
                    integral_orbitals(3),integral_orbitals(2),2,exchange_breit_coefficients(9,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(3),integral_orbitals(2), &
                    integral_orbitals(1),integral_orbitals(4),2,exchange_breit_coefficients(10,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(4),integral_orbitals(1), &
                    integral_orbitals(2),integral_orbitals(3),2,exchange_breit_coefficients(11,breit_rank_index))
                call talk(bra_csf_index,ket_csf_index,next_breit_radial_rank,integral_orbitals(2),integral_orbitals(3), &
                    integral_orbitals(4),integral_orbitals(1),2,exchange_breit_coefficients(12,breit_rank_index))
            end do
        end if
        return

    contains
        subroutine report_rank_limit()
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN EL53  PMGG RAGG KRA1=",i100)')rank_table_position
            stop
        end subroutine report_rank_limit
    end subroutine el53
end module four_shell_transfer_14_23
