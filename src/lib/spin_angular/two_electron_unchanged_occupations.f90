! Evaluate Coulomb or Breit terms with unchanged subshell occupations.
module two_electron_unchanged_occupations
    implicit none (type, external)

    ! Module interface
    private
    public :: accumulate_two_electron_unchanged_occupations
contains
    subroutine accumulate_two_electron_unchanged_occupations(bra_csf_index,ket_csf_index,shell_position_a,shell_position_b, &
        & check_recoupling,interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, half, angular_tolerance
        use csf_pair_angular_state,             only: peel_orbital_indices, peel_orbital_count
        use orbital_csf_basis,           only: orbital_kappa
        use jj_tensor_coupling_tables, only: ket_shell_values_1, bra_shell_values_1, ket_shell_codes_1, ket_shell_codes_2, &
            bra_shell_codes_1
        use scalar_interaction_recoupling_selection, only: select_scalar_interaction_recoupling
        use two_subshell_recoupling, only: recouple_two_subshells
        use jj_active_subshell_states, only: prepare_jj_active_subshell_states
        use angular_momentum_triangle, only: satisfies_angular_triangle
        use shared_tensor_rank_range, only: shared_tensor_rank_bounds
        use wigner_6j_triangle_selection, only: allows_wigner_6j_triangles
        use breit_radial_rank_ranges, only: breit_rank_bounds
        use mcp_coefficient_output, only: buffer_mcp_coefficient
        use coulomb_angular_factor, only: coulomb_angular_factor_value
        use jj_ww_scalar_matrix_element, only: evaluate_jj_scalar_ww
        use wigner_6j, only: wigner_6j_value
        use breit_radial_integral_coefficients, only: calculate_breit_radial_integral_coefficients
        use breit_angular_coefficient_buffer, only: collect_breit_angular_coefficient
        use subshell_product_w_w, only: evaluate_subshell_product_w_w
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,shell_position_a,shell_position_b,check_recoupling, &
            interaction_type

        ! Local declarations
        integer :: work_index,orbital_a,orbital_b,coupling_allowed,first_direct_rank_weight, &
            first_interaction_rank_weight,first_saved_rank_weight,last_direct_rank_weight, &
            last_interaction_rank_weight,last_saved_rank_weight,two_rank_range_size,two_direct_rank_plus_one, &
            two_multipole_rank_plus_one,two_coupled_rank_plus_one,two_precomputed_rank_plus_one,&
            phase_sum,coupled_rank,direct_integral_type,exchange_integral_type,multipole_rank, &
            rank_table_position,abs_kappa_1,abs_kappa_2,breit_component,breit_rank_index, &
            breit_radial_rank,first_direct_nu,direct_nu_count,   &
            first_exchange_nu,exchange_nu_count,next_breit_radial_rank
        integer :: two_orbital_j(2)
        integer :: integral_orbitals(4),two_kappas(4),two_j_plus_one(4)
        real(real64)          :: q_projection_1,q_projection_2,q_projection_3,q_projection_4,coupled_term, &
            angular_sum,coulomb_angular_factor,angular_coefficient,six_j_value,recoupling_factor,subshell_matrix_element
        real(real64)   :: breit_factors(12)
        real(real64)   :: recoupling_by_rank(30),subshell_element_by_rank(30)
        real(real64) :: exchange_breit_coefficients(7, 20)
        if (shell_position_a == shell_position_b) then
            if(check_recoupling /= 0) then
                call select_scalar_interaction_recoupling(shell_position_a,shell_position_a,shell_position_a,shell_position_a,0, &
                    & coupling_allowed)
                if(coupling_allowed == 0)return
            end if
            call prepare_jj_active_subshell_states(shell_position_a,shell_position_a,shell_position_a,shell_position_a,1)
            q_projection_1=half
            q_projection_2=half
            q_projection_3=-half
            q_projection_4=-half
            orbital_a=peel_orbital_indices(shell_position_a)
            two_orbital_j(1)=ket_shell_codes_1(3)
            first_interaction_rank_weight=shared_tensor_rank_bounds(two_orbital_j(1),two_orbital_j(1),two_orbital_j(1), &
                & two_orbital_j(1), &
                two_rank_range_size)+1
            if(two_rank_range_size <= 0) return
            last_interaction_rank_weight=first_interaction_rank_weight+two_rank_range_size-1
            abs_kappa_1=(two_orbital_j(1)+1)/2
            first_direct_rank_weight=first_interaction_rank_weight
            last_direct_rank_weight=last_interaction_rank_weight
            if (interaction_type == 2) then
                integral_orbitals(1)=orbital_a
                integral_orbitals(2)=orbital_a
                integral_orbitals(3)=orbital_a
                integral_orbitals(4)=orbital_a
                two_kappas(1)=2*orbital_kappa(integral_orbitals(1))
                two_kappas(2)=2*orbital_kappa(integral_orbitals(2))
                two_kappas(3)=2*orbital_kappa(integral_orbitals(3))
                two_kappas(4)=2*orbital_kappa(integral_orbitals(4))
                two_j_plus_one(1)=iabs(two_kappas(1))
                two_j_plus_one(2)=iabs(two_kappas(2))
                two_j_plus_one(3)=iabs(two_kappas(3))
                two_j_plus_one(4)=iabs(two_kappas(4))
                call breit_rank_bounds(integral_orbitals,two_kappas,two_j_plus_one,first_direct_nu,direct_nu_count, &
                    & first_exchange_nu, &
                    exchange_nu_count,direct_integral_type,exchange_integral_type)
                if(direct_integral_type <= 0)return
            end if
            do two_multipole_rank_plus_one=first_direct_rank_weight,last_direct_rank_weight,2
                multipole_rank=(two_multipole_rank_plus_one-1)/2
                if (interaction_type == 1) then
                    call coulomb_angular_factor_value(abs_kappa_1,abs_kappa_1,abs_kappa_1,abs_kappa_1,ket_shell_codes_1(5), &
                        & ket_shell_codes_1(5), &
                        ket_shell_codes_1(5),ket_shell_codes_1(5),multipole_rank,coulomb_angular_factor)
                    if(dabs(coulomb_angular_factor) < angular_tolerance)cycle
                    coulomb_angular_factor=-coulomb_angular_factor*half
                end if

                ! Sum the 6j-transformed subshell factors over intermediate tensor ranks.
                angular_sum=zero
                do two_coupled_rank_plus_one=first_interaction_rank_weight,last_interaction_rank_weight,2
                    coupled_rank=(two_coupled_rank_plus_one-1)/2
                    if(allows_wigner_6j_triangles(two_orbital_j(1),two_orbital_j(1),multipole_rank*2,two_orbital_j(1), &
                        & two_orbital_j(1), &
                        coupled_rank*2) == 0)cycle
                    call evaluate_jj_scalar_ww(bra_shell_codes_1,bra_shell_values_1,ket_shell_codes_1,ket_shell_values_1, &
                        & coupled_rank, &
                        q_projection_1,q_projection_2,q_projection_3,q_projection_4,coupled_term)
                    if(dabs(coupled_term) < angular_tolerance)cycle
                    call wigner_6j_value(two_orbital_j(1),two_orbital_j(1),multipole_rank*2,two_orbital_j(1),two_orbital_j(1), &
                        coupled_rank*2,0,six_j_value)
                    coupled_term=coupled_term*six_j_value*dsqrt(dble(two_coupled_rank_plus_one))
                    phase_sum=bra_shell_codes_1(3)+coupled_rank+multipole_rank
                    if((phase_sum/2)*2 /= phase_sum)coupled_term=-coupled_term
                    angular_sum=angular_sum+coupled_term
                end do
                if (interaction_type == 1) then
                    angular_coefficient=angular_sum*coulomb_angular_factor
                    angular_coefficient=angular_coefficient/dsqrt(dble(bra_shell_codes_1(6)+1))
                    if(dabs(angular_coefficient) > angular_tolerance)call buffer_mcp_coefficient(bra_csf_index,ket_csf_index, &
                        & orbital_a, &
                        orbital_a,orbital_a,orbital_a,multipole_rank,angular_coefficient)
                else if (interaction_type == 2) then
                    breit_rank_index=(multipole_rank-first_direct_nu)/2+1
                    if(((multipole_rank-first_direct_nu)/2)*2 == (multipole_rank-first_direct_nu)) then
                        call calculate_breit_radial_integral_coefficients(breit_factors, integral_orbitals, &
                            two_kappas, multipole_rank, multipole_rank, 3, 1)
                        if(dabs(breit_factors(1)) > angular_tolerance) then
                            angular_coefficient =-half*breit_factors(1)*angular_sum/dsqrt(dble(bra_shell_codes_1(6)+1))
                            if(dabs(angular_coefficient) > angular_tolerance)call collect_breit_angular_coefficient( &
                                & bra_csf_index,ket_csf_index, &
                                multipole_rank,orbital_a,orbital_a,orbital_a,orbital_a,4,angular_coefficient)
                        end if
                    end if
                end if
            end do
            return
        end if
        if(peel_orbital_count <= 1)return
        if(check_recoupling /= 0) then
            call select_scalar_interaction_recoupling(shell_position_a,shell_position_b,shell_position_b,shell_position_b,1, &
                & coupling_allowed)
            if(coupling_allowed == 0)return
        end if
        orbital_a=peel_orbital_indices(shell_position_a)
        orbital_b=peel_orbital_indices(shell_position_b)
        q_projection_1=half
        q_projection_2=-half
        q_projection_3=half
        q_projection_4=-half
        call prepare_jj_active_subshell_states(shell_position_a,shell_position_b,shell_position_a,shell_position_a,2)
        two_orbital_j(1)=ket_shell_codes_1(3)
        two_orbital_j(2)=ket_shell_codes_2(3)
        abs_kappa_1=(two_orbital_j(1)+1)/2
        abs_kappa_2=(two_orbital_j(2)+1)/2
        first_direct_rank_weight=shared_tensor_rank_bounds(two_orbital_j(1),two_orbital_j(1),two_orbital_j(2),two_orbital_j(2), &
            two_rank_range_size)+1
        if(two_rank_range_size <= 0)return
        last_direct_rank_weight=first_direct_rank_weight+two_rank_range_size-1
        first_saved_rank_weight=first_direct_rank_weight
        last_saved_rank_weight=last_direct_rank_weight
        do two_precomputed_rank_plus_one=first_direct_rank_weight,last_direct_rank_weight,2
            multipole_rank=(two_precomputed_rank_plus_one-1)/2
            rank_table_position=multipole_rank+1
            if(rank_table_position > 30) call report_rank_limit()
            subshell_element_by_rank(rank_table_position)=zero

            ! Cache the subshell recoupling factor at each intermediate rank for later orbital pairings.
            recoupling_by_rank(rank_table_position)=zero
            call recouple_two_subshells(shell_position_a,shell_position_b,multipole_rank*2,0,coupling_allowed,recoupling_factor)
            if(coupling_allowed == 0) cycle
            call evaluate_subshell_product_w_w(multipole_rank,multipole_rank,q_projection_1,q_projection_2,q_projection_3, &
                & q_projection_4, &
                subshell_matrix_element)
            if(dabs(subshell_matrix_element) < angular_tolerance) cycle
            subshell_element_by_rank(rank_table_position)=subshell_matrix_element
            call recouple_two_subshells(shell_position_a,shell_position_b,multipole_rank*2,1,coupling_allowed,recoupling_factor)
            recoupling_by_rank(rank_table_position)=recoupling_factor
        end do
        if (interaction_type == 2) then
            integral_orbitals(1)=orbital_a
            integral_orbitals(2)=orbital_b
            integral_orbitals(3)=orbital_a
            integral_orbitals(4)=orbital_b
            two_kappas(1)=2*orbital_kappa(integral_orbitals(1))
            two_kappas(2)=2*orbital_kappa(integral_orbitals(2))
            two_kappas(3)=2*orbital_kappa(integral_orbitals(3))
            two_kappas(4)=2*orbital_kappa(integral_orbitals(4))
            two_j_plus_one(1)=iabs(two_kappas(1))
            two_j_plus_one(2)=iabs(two_kappas(2))
            two_j_plus_one(3)=iabs(two_kappas(3))
            two_j_plus_one(4)=iabs(two_kappas(4))
            call breit_rank_bounds(integral_orbitals,two_kappas,two_j_plus_one,first_direct_nu,direct_nu_count,first_exchange_nu, &
                exchange_nu_count,direct_integral_type,exchange_integral_type)
            do work_index=1,20
                exchange_breit_coefficients(1,work_index)=zero
                exchange_breit_coefficients(2,work_index)=zero
                exchange_breit_coefficients(3,work_index)=zero
                exchange_breit_coefficients(4,work_index)=zero
                exchange_breit_coefficients(5,work_index)=zero
                exchange_breit_coefficients(6,work_index)=zero
                exchange_breit_coefficients(7,work_index)=zero
            end do
            if(direct_integral_type == 0 .and. exchange_integral_type == 0)return
        end if
        do two_direct_rank_plus_one=first_direct_rank_weight,last_direct_rank_weight,2
            multipole_rank=(two_direct_rank_plus_one-1)/2
            rank_table_position=multipole_rank+1
            if(rank_table_position > 30) call report_rank_limit()
            if (interaction_type == 1) then
                call coulomb_angular_factor_value(abs_kappa_1,abs_kappa_2,abs_kappa_1,abs_kappa_2,ket_shell_codes_1(5), &
                    & ket_shell_codes_2(5), &
                    ket_shell_codes_1(5),ket_shell_codes_2(5),multipole_rank,coupled_term)
                if(dabs(coupled_term) < angular_tolerance) cycle
                coupled_term=coupled_term*recoupling_by_rank(rank_table_position)
                if(dabs(coupled_term) < angular_tolerance) cycle
                coupled_term=coupled_term*subshell_element_by_rank(rank_table_position)
                if(dabs(coupled_term) < angular_tolerance) cycle
                coupled_term=coupled_term/dsqrt(dble(two_direct_rank_plus_one))
                if(dabs(coupled_term) > angular_tolerance) call buffer_mcp_coefficient(bra_csf_index,ket_csf_index,orbital_a, &
                    & orbital_b, &
                    orbital_a,orbital_b,multipole_rank,coupled_term)
            else if (interaction_type == 2) then
                breit_rank_index=(multipole_rank-first_direct_nu)/2+1
                if(((multipole_rank-first_direct_nu)/2)*2 == (multipole_rank-first_direct_nu)) then
                    call calculate_breit_radial_integral_coefficients(breit_factors, integral_orbitals, &
                        two_kappas, multipole_rank, multipole_rank, 3, 1)
                    if(dabs(breit_factors(1)) > angular_tolerance) then
                        angular_coefficient=breit_factors(1)*recoupling_by_rank(rank_table_position)* &
                            subshell_element_by_rank(rank_table_position)/dsqrt(dble(two_direct_rank_plus_one))
                        if(dabs(angular_coefficient) > angular_tolerance)call collect_breit_angular_coefficient(bra_csf_index, &
                            & ket_csf_index, &
                            multipole_rank,orbital_a,orbital_a,orbital_b,orbital_b,4,angular_coefficient)
                    end if
                end if
            end if
        end do
        first_interaction_rank_weight=shared_tensor_rank_bounds(two_orbital_j(1),two_orbital_j(2),two_orbital_j(1), &
            & two_orbital_j(2), &
            two_rank_range_size)+1
        if(two_rank_range_size <= 0) return
        last_interaction_rank_weight=first_interaction_rank_weight+two_rank_range_size-1
        do two_multipole_rank_plus_one=first_interaction_rank_weight,last_interaction_rank_weight,2
            multipole_rank=(two_multipole_rank_plus_one-1)/2
            if(multipole_rank > 30) call report_rank_limit()
            if (interaction_type == 1) then
                call coulomb_angular_factor_value(abs_kappa_1,abs_kappa_2,abs_kappa_2,abs_kappa_1,ket_shell_codes_1(5), &
                    & ket_shell_codes_2(5), &
                    ket_shell_codes_2(5),ket_shell_codes_1(5),multipole_rank,coulomb_angular_factor)
                if(dabs(coulomb_angular_factor) < angular_tolerance) cycle
            end if
            angular_sum=zero
            do two_coupled_rank_plus_one=first_saved_rank_weight,last_saved_rank_weight,2
                coupled_rank=(two_coupled_rank_plus_one-1)/2
                rank_table_position=coupled_rank+1
                if(rank_table_position > 30) call report_rank_limit()
                coupled_term=recoupling_by_rank(rank_table_position)
                if(dabs(coupled_term) < angular_tolerance) cycle
                coupled_term=coupled_term*subshell_element_by_rank(rank_table_position)
                if(dabs(coupled_term) < angular_tolerance) cycle
                if(allows_wigner_6j_triangles(two_orbital_j(1),two_orbital_j(2),multipole_rank*2,two_orbital_j(2), &
                    & two_orbital_j(1), &
                    coupled_rank*2) == 0)cycle
                call wigner_6j_value(two_orbital_j(1),two_orbital_j(2),multipole_rank*2,two_orbital_j(2),two_orbital_j(1), &
                    coupled_rank*2,0,six_j_value)
                coupled_term=coupled_term*six_j_value*dsqrt(dble(two_coupled_rank_plus_one))
                angular_sum=angular_sum+coupled_term
            end do
            if (interaction_type == 1) then
                angular_coefficient=coulomb_angular_factor*angular_sum
                if(dabs(angular_coefficient) > angular_tolerance)call buffer_mcp_coefficient(bra_csf_index,ket_csf_index, &
                    & orbital_a, &
                    orbital_b,orbital_b,orbital_a,multipole_rank,angular_coefficient)
            else if (interaction_type == 2) then
                breit_radial_rank=multipole_rank
                if(((breit_radial_rank-first_exchange_nu)/2)*2 == (breit_radial_rank-first_exchange_nu)) then
                    if((satisfies_angular_triangle(two_j_plus_one(1),two_j_plus_one(4), &
                        & breit_radial_rank+breit_radial_rank+1) /= 0) .and.  &
                        (satisfies_angular_triangle(two_j_plus_one(2),two_j_plus_one(3), &
                            & breit_radial_rank+breit_radial_rank+1) /= 0)) then
                        if(breit_radial_rank > 0) then
                            breit_rank_index=(breit_radial_rank-first_exchange_nu)/2+1
                            call calculate_breit_radial_integral_coefficients(breit_factors, integral_orbitals, &
                                two_kappas, breit_radial_rank, multipole_rank, 4, 2)
                            do breit_component = 1,3
                                exchange_breit_coefficients(breit_component,breit_rank_index)= &
                                    exchange_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
                                    breit_factors(breit_component)
                            end do
                        end if
                    end if
                end if
                breit_radial_rank=multipole_rank+1
                if(((breit_radial_rank-first_exchange_nu)/2)*2 == (breit_radial_rank-first_exchange_nu)) then
                    if((satisfies_angular_triangle(two_j_plus_one(1),two_j_plus_one(4), &
                        & breit_radial_rank+breit_radial_rank-1) /= 0) .and.  &
                        (satisfies_angular_triangle(two_j_plus_one(2),two_j_plus_one(3), &
                            & breit_radial_rank+breit_radial_rank-1) /= 0)) then
                        if(breit_radial_rank >= 0) then
                            breit_rank_index=(breit_radial_rank-first_exchange_nu)/2+1
                            if(breit_rank_index <= exchange_nu_count) then
                                call calculate_breit_radial_integral_coefficients(breit_factors, integral_orbitals, &
                                    two_kappas, breit_radial_rank, multipole_rank, 4, 2)
                                do breit_component = 1,3
                                    exchange_breit_coefficients(breit_component,breit_rank_index)= &
                                        exchange_breit_coefficients(breit_component,breit_rank_index)+angular_sum* &
                                        breit_factors(breit_component)
                                end do
                            end if
                        end if
                    end if
                end if
                breit_radial_rank=multipole_rank-1
                if(((breit_radial_rank-first_exchange_nu)/2)*2 == (breit_radial_rank-first_exchange_nu)) then
                    if((satisfies_angular_triangle(two_j_plus_one(1),two_j_plus_one(4), &
                        & breit_radial_rank+breit_radial_rank+3) /= 0) .and.  &
                        (satisfies_angular_triangle(two_j_plus_one(2),two_j_plus_one(3), &
                            & breit_radial_rank+breit_radial_rank+3) /= 0)) then
                        if(breit_radial_rank >= 0) then
                            breit_rank_index=(breit_radial_rank-first_exchange_nu)/2+1
                            if(breit_rank_index < exchange_nu_count) then
                                call calculate_breit_radial_integral_coefficients(breit_factors, integral_orbitals, &
                                    two_kappas, breit_radial_rank, multipole_rank, 4, 2)
                                do breit_component = 1,7
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
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,breit_radial_rank,orbital_b,orbital_a, &
                    & orbital_b,orbital_a,5, &
                    exchange_breit_coefficients(1,breit_rank_index))
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,breit_radial_rank,orbital_a,orbital_b, &
                    & orbital_b,orbital_a,5, &
                    exchange_breit_coefficients(2,breit_rank_index))
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,breit_radial_rank,orbital_a,orbital_b, &
                    & orbital_a,orbital_b,5, &
                    exchange_breit_coefficients(3,breit_rank_index))
                if(breit_rank_index == exchange_nu_count) cycle
                next_breit_radial_rank=breit_radial_rank+1
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,next_breit_radial_rank,orbital_a,orbital_b, &
                    & orbital_a,orbital_b,6, &
                    exchange_breit_coefficients(4,breit_rank_index))
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,next_breit_radial_rank,orbital_b,orbital_a, &
                    & orbital_b,orbital_a,6, &
                    exchange_breit_coefficients(5,breit_rank_index))
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,next_breit_radial_rank,orbital_a,orbital_b, &
                    & orbital_b,orbital_a,6, &
                    exchange_breit_coefficients(6,breit_rank_index))
                call collect_breit_angular_coefficient(bra_csf_index,ket_csf_index,next_breit_radial_rank,orbital_b,orbital_a, &
                    & orbital_a,orbital_b,6, &
                    exchange_breit_coefficients(7,breit_rank_index))
            end do
        end if
        return

    contains
        subroutine report_rank_limit()
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN EL1  PMGG RAGG")')
            stop
        end subroutine report_rank_limit
    end subroutine accumulate_two_electron_unchanged_occupations
end module two_electron_unchanged_occupations
