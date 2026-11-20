! Dispatch Coulomb or Breit coefficients by symmetry and CSF occupation changes.
module two_electron_angular_coefficients
    implicit none (type, external)

    ! Module interface
    private
    public :: rkco_gg
contains
    subroutine rkco_gg (bra_csf_index,ket_csf_index,coefficient_consumer,include_closed_shell_terms, &
        interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use csf_pair_angular_state, only: bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2, &
            bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2, peel_orbital_indices => jlist, &
            core_orbital_indices => klist, max_orbitals, core_orbital_count => ncore, peel_orbital_count => npeel, &
            bra_occupations => nq1, ket_occupations => nq2
        use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa
        use angular_recoupling_backup, only: saved_bra_j_weights => jc1s, saved_ket_j_weights => jc2s, &
            saved_peel_orbitals => jlis, max_orbitals
        use csf_quantum_queries, only: itjpo, ispar
        use two_electron_csf_filter, only: ichkq2
        use csf_pair_quantum_numbers, only: setqna
        use two_electron_unchanged_occupations, only: el1
        use two_electron_pair_transfer, only: el2
        use two_electron_single_transfer, only: el3
        use two_electron_three_shell_transfer, only: el4
        use two_electron_four_shell_transfer, only: el5
        use closed_shell_mcp, only: cord_callback => cord
        implicit none (type, external)

        ! Arguments
        procedure(cord_callback) :: coefficient_consumer
        integer, intent(in) :: bra_csf_index,ket_csf_index,include_closed_shell_terms,interaction_type

        ! Local declarations
        integer :: shell_position,bra_orbital,ket_orbital,total_occupation_difference,occupation_pattern_code, &
            shift_step,first_orbital,orbital_index,core_orbital,shift_count,orbital_list_selector,first_bra_orbital, &
            first_ket_orbital, &
            bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,peel_orbital,peel_position,saved_bra_shell_1, &
            saved_ket_shell_1,swap_shell_position,occupation_difference, &
            last_outer_position,partner_shell_position,outer_peel_position,partner_orbital, &
            second_shell_position,first_shell_position,peel_coupling_count

        ! Coulomb and Breit are even scalars, so the two CSFs must have equal J and parity.
        if (itjpo (bra_csf_index) /= itjpo (ket_csf_index)) return
        if (ispar (bra_csf_index) /= ispar (ket_csf_index)) return
        if (ichkq2(bra_csf_index,ket_csf_index) == 0) return
        call setqna (bra_csf_index,ket_csf_index)
        if (orbital_count < 1) then
            print *, "RKCO_GG: No subshells."
            stop
        end if
        if (peel_orbital_count /= 0 .and. peel_orbital_count /= 1) then
            total_occupation_difference = 0
            bra_shell_1 = 0
            bra_shell_2 = 0
            ket_shell_1 = 0
            ket_shell_2 = 0

            ! Record both the total absolute occupation difference and the pattern of one- and two-electron changes.
            ! A changed shell contributes 1 for a one-electron change or 20 for a two-electron change.
            ! With two electrons transferred, codes 4, 22 and 40 identify four, three and two shells.
            occupation_pattern_code = 0
            do peel_position = 1,peel_orbital_count
                peel_orbital = peel_orbital_indices(peel_position)
                occupation_difference = bra_occupations(peel_orbital) - ket_occupations(peel_orbital)
                if (iabs (occupation_difference) > 2) return
                if (occupation_difference < 0) then
                    if (occupation_difference+1 > 0) then
                        cycle
                    else if (occupation_difference+1 == 0) then
                        if (ket_shell_1 > 0) then
                            ket_shell_2 = peel_position
                        else
                            ket_shell_1 = peel_position
                        end if
                        total_occupation_difference = total_occupation_difference+1
                        occupation_pattern_code=1+occupation_pattern_code
                    else if (occupation_difference+1 < 0) then
                        ket_shell_1 = peel_position
                        total_occupation_difference = total_occupation_difference+2
                        occupation_pattern_code=20+occupation_pattern_code
                    end if
                else
                    if (occupation_difference-1 > 0) then
                        bra_shell_1 = peel_position
                        total_occupation_difference = total_occupation_difference+2
                        occupation_pattern_code=20+occupation_pattern_code
                        cycle
                    else if (occupation_difference-1 == 0) then
                        if (bra_shell_1 > 0) then
                            bra_shell_2 = peel_position
                        else
                            bra_shell_1 = peel_position
                        end if
                        total_occupation_difference = total_occupation_difference+1
                        occupation_pattern_code=1+occupation_pattern_code
                        cycle
                    else
                        cycle
                    end if
                end if
            end do
            if (total_occupation_difference > 4) return

            ! Two transferred electrons can occupy two, three or four distinct subshells.
            if (total_occupation_difference == 4) then
                if (bra_shell_2 == 0) then
                    bra_shell_2 = bra_shell_1
                end if
                if (ket_shell_2 == 0) then
                    ket_shell_2 = ket_shell_1
                end if
                if(occupation_pattern_code/=40) then
                    if(occupation_pattern_code/=22) then
                        call el5(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,interaction_type)
                    else
                        call el4(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,interaction_type)
                    end if
                else
                    call el2(bra_csf_index,ket_csf_index,bra_shell_1,ket_shell_1,interaction_type)
                end if
                return
            end if
            select case (total_occupation_difference)
            case (0)
                ! Equal CSFs proceed to the diagonal terms below.
                last_outer_position = 0
                if (bra_csf_index /= ket_csf_index) last_outer_position = peel_orbital_count
            case (2)
                ! The changed shell pair is already known. Loop only over occupied partners.
                last_outer_position = 1
            case default
                write (*,'("RKCO_GG: Error.")')
                stop
            end select
            do outer_peel_position = 1,last_outer_position
                if (total_occupation_difference /= 2) then
                    bra_shell_1 = outer_peel_position
                    ket_shell_1 = outer_peel_position
                end if
                saved_bra_shell_1 = bra_shell_1
                saved_ket_shell_1 = ket_shell_1
                first_bra_orbital = peel_orbital_indices(bra_shell_1)
                first_ket_orbital = peel_orbital_indices(ket_shell_1)
                do partner_shell_position = outer_peel_position,peel_orbital_count
                    partner_orbital = peel_orbital_indices(partner_shell_position)
                    if (bra_occupations(partner_orbital)*ket_occupations(partner_orbital) == 0) cycle
                    bra_shell_2 = partner_shell_position
                    ket_shell_2 = partner_shell_position
                    bra_shell_1 = saved_bra_shell_1
                    ket_shell_1 = saved_ket_shell_1
                    if (bra_shell_1-bra_shell_2 > 0) then
                        swap_shell_position = bra_shell_2
                        bra_shell_2 = bra_shell_1
                        bra_shell_1 = swap_shell_position
                    else if (bra_shell_1-bra_shell_2 == 0) then
                        bra_orbital = peel_orbital_indices(bra_shell_2)
                        if (bra_occupations(bra_orbital) <= 1) cycle
                    end if
                    if (ket_shell_1-ket_shell_2 > 0) then
                        swap_shell_position = ket_shell_2
                        ket_shell_2 = ket_shell_1
                        ket_shell_1 = swap_shell_position
                    else if (ket_shell_1-ket_shell_2 == 0) then
                        ket_orbital = peel_orbital_indices(ket_shell_2)
                        if (ket_occupations(ket_orbital) <= 1) cycle
                    end if
                    if(total_occupation_difference/=0) then
                        call el3(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,interaction_type)
                    else
                        call el1(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,1,interaction_type)
                    end if
                end do
                if ((total_occupation_difference == 0) .and. (core_orbital_count == 0)) cycle
                if ((core_orbital_count == 0) .or. (orbital_kappa(first_bra_orbital) /= orbital_kappa(first_ket_orbital))) &
                    return
                peel_coupling_count = peel_orbital_count-1
                do shell_position = 1,peel_orbital_count
                    saved_peel_orbitals(shell_position) = peel_orbital_indices(shell_position)
                end do
                do shell_position = 1,peel_coupling_count
                    saved_bra_j_weights(shell_position) = bra_intermediate_couplings(shell_position)
                    saved_ket_j_weights(shell_position) = ket_intermediate_couplings(shell_position)
                end do

                ! Insert a core orbital as a partner and restore the saved peel coupling after evaluation.
                do partner_shell_position = 1,core_orbital_count
                    core_orbital = core_orbital_indices(partner_shell_position)
                    do shell_position = 1,peel_orbital_count
                        orbital_index = peel_orbital_indices(shell_position)
                        if (core_orbital < orbital_index) exit
                    end do
                    shift_count = peel_orbital_count-shell_position+1
                    do shift_step = 1,shift_count
                        peel_orbital_indices(peel_orbital_count+2-shift_step) = peel_orbital_indices(peel_orbital_count+1- &
                            shift_step)
                        if (peel_orbital_count == shift_step) exit
                        bra_intermediate_couplings(peel_orbital_count+1-shift_step) = &
                            bra_intermediate_couplings(peel_orbital_count-shift_step)
                        ket_intermediate_couplings(peel_orbital_count+1-shift_step) = &
                            ket_intermediate_couplings(peel_orbital_count-shift_step)
                    end do
                    if (shell_position < 3) then
                        first_orbital = peel_orbital_indices(1)
                        bra_intermediate_couplings(1) = bra_subshell_quantum_codes(3,first_orbital)
                        ket_intermediate_couplings(1) = ket_subshell_quantum_codes(3,first_orbital)
                    else
                        bra_intermediate_couplings(shell_position-1) = bra_intermediate_couplings(shell_position-2)
                        ket_intermediate_couplings(shell_position-1) = ket_intermediate_couplings(shell_position-2)
                    end if
                    peel_orbital_indices(shell_position) = core_orbital
                    bra_shell_1 = saved_bra_shell_1
                    if (saved_bra_shell_1 >= shell_position) bra_shell_1 = bra_shell_1+1
                    bra_shell_2 = shell_position
                    ket_shell_1 = saved_ket_shell_1
                    if (saved_ket_shell_1 >= shell_position) ket_shell_1 = ket_shell_1+1
                    ket_shell_2 = shell_position
                    if (bra_shell_1-bra_shell_2 > 0) then
                        swap_shell_position = bra_shell_2
                        bra_shell_2 = bra_shell_1
                        bra_shell_1 = swap_shell_position
                    end if
                    if (ket_shell_1-ket_shell_2 > 0) then
                        swap_shell_position = ket_shell_2
                        ket_shell_2 = ket_shell_1
                        ket_shell_1 = swap_shell_position
                    end if
                    peel_orbital_count = peel_orbital_count+1
                    if(total_occupation_difference/=0) then
                        if(occupation_pattern_code/=40) then
                            if(occupation_pattern_code/=2) then
                                call ensure_debug_output()
                                write(angular_debug_unit,'("   rie zymes 38 atv N=N-N  !!!!!!")')
                                return
                            else
                                call el3(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
                                    interaction_type)
                            end if
                        else
                            call ensure_debug_output()
                            write(angular_debug_unit,'("   rie zymes 38?? atv N=N-N  !!!!!!")')
                            return
                        end if
                    else
                        call el1(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,1,interaction_type)
                    end if
                    peel_orbital_count = peel_orbital_count-1
                    peel_coupling_count = peel_orbital_count-1
                    do shell_position = 1,peel_orbital_count
                        peel_orbital_indices(shell_position) = saved_peel_orbitals(shell_position)
                    end do
                    do shell_position = 1,peel_coupling_count
                        bra_intermediate_couplings(shell_position)  = saved_bra_j_weights(shell_position)
                        ket_intermediate_couplings(shell_position)  = saved_ket_j_weights(shell_position)
                    end do
                end do
            end do
            if (last_outer_position /= 0) return
        end if

        do second_shell_position = 1,peel_orbital_count
            partner_orbital = peel_orbital_indices(second_shell_position)
            bra_shell_2 = second_shell_position
            ket_shell_2 = second_shell_position
            do first_shell_position = 1,second_shell_position
                bra_shell_1 = first_shell_position
                if (bra_shell_1 == bra_shell_2) then
                    if (bra_occupations(partner_orbital) <= 1) cycle
                end if
                ket_shell_1 = bra_shell_1
                if(bra_csf_index /= ket_csf_index) then
                    if(peel_orbital_count == 1 .and. bra_occupations(partner_orbital) == ket_occupations(partner_orbital) &
                        ) then
                        call el1(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,1,interaction_type)
                    else if(occupation_pattern_code/=2) then
                        call ensure_debug_output()
                        write(angular_debug_unit,'("   rie zymes 45 atv N=N-N  !!!!!!")')
                        call ensure_debug_output()
                        write(angular_debug_unit,*)"JA,JB,JA1,JB1",bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2
                        call ensure_debug_output()
                        write(angular_debug_unit,*)"IDQG IDQ",occupation_pattern_code,total_occupation_difference
                        return
                    else
                        call el3(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
                            interaction_type)
                    end if
                else
                    call el1(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,0,interaction_type)
                end if
            end do
        end do

        ! On the CSF diagonal, use the callback for contributions with closed core subshells.
        if (include_closed_shell_terms < 1) return
        if (core_orbital_count == 0) return
        if(bra_csf_index/=ket_csf_index) return
        do second_shell_position = 1,core_orbital_count
            bra_shell_2 = second_shell_position
            ket_shell_2 = second_shell_position
            orbital_list_selector = 2
            do first_shell_position = 1,second_shell_position
                bra_shell_1 = first_shell_position
                ket_shell_1 = first_shell_position
                call coefficient_consumer (bra_csf_index,ket_csf_index,bra_shell_1,orbital_list_selector,bra_shell_2)
            end do
            if (peel_orbital_count == 0) cycle
            orbital_list_selector = 1
            do first_shell_position = 1,peel_orbital_count
                bra_shell_1 = first_shell_position
                ket_shell_1 = first_shell_position
                call coefficient_consumer (bra_csf_index,ket_csf_index,bra_shell_1,orbital_list_selector,bra_shell_2)
            end do
        end do
        return
    end subroutine rkco_gg
end module two_electron_angular_coefficients
