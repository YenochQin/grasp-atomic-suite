! Dispatch two-electron terms with a net transfer of one electron.
module two_electron_single_transfer
    implicit none (type, external)

    ! Module interface
    private
    public :: accumulate_two_electron_single_transfer
contains
    subroutine accumulate_two_electron_single_transfer(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1, &
        & ket_shell_2,interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use csf_pair_angular_state,         only:  peel_orbital_count
        use single_transfer_triple_on_first_shell, only: accumulate_single_transfer_triple_on_first_shell
        use single_transfer_triple_on_second_shell, only: accumulate_single_transfer_triple_on_second_shell
        use single_transfer_with_spectator_shell, only: accumulate_single_transfer_with_spectator_shell
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
            interaction_type
        if(peel_orbital_count <= 1)return
        if(bra_shell_2 == ket_shell_2) then
            if(bra_shell_1 == bra_shell_2.or.ket_shell_1 == bra_shell_2) then
                if(bra_shell_1 == ket_shell_1) call report_invalid_shell_transfer()
                if(ket_shell_1 /= bra_shell_2) then
                    call accumulate_single_transfer_triple_on_second_shell(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                else
                    call accumulate_single_transfer_triple_on_first_shell(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                end if
            else
                call accumulate_single_transfer_with_spectator_shell(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1, &
                    & bra_shell_2,1,bra_shell_1,bra_shell_2, &
                    ket_shell_1,ket_shell_2,interaction_type)
            end if
            return
        else if(bra_shell_1 == ket_shell_1) then
            if(bra_shell_2 == bra_shell_1.or.ket_shell_2 == bra_shell_1) then
                if(bra_shell_2 == ket_shell_2) call report_invalid_shell_transfer()
                if(ket_shell_2 /= bra_shell_1) then
                    call accumulate_single_transfer_triple_on_second_shell(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_2, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                else
                    call accumulate_single_transfer_triple_on_first_shell(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_2, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                end if
            else
                call accumulate_single_transfer_with_spectator_shell(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_2, &
                    & bra_shell_1,1,bra_shell_1,bra_shell_2, &
                    ket_shell_1,ket_shell_2,interaction_type)
            end if
            return
        else if(bra_shell_1 == ket_shell_2) then
            if(bra_shell_2 == bra_shell_1.or.ket_shell_1 == bra_shell_1) then
                if(bra_shell_2 == ket_shell_1) call report_invalid_shell_transfer()
                if(ket_shell_1 /= ket_shell_2) then
                    call accumulate_single_transfer_triple_on_second_shell(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_2, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                else
                    call accumulate_single_transfer_triple_on_first_shell(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_2, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                end if
            else
                call accumulate_single_transfer_with_spectator_shell(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_2, &
                    & bra_shell_1,2,bra_shell_1,bra_shell_2, &
                    ket_shell_2,ket_shell_1,interaction_type)
            end if
            return
        else if(bra_shell_2 == ket_shell_1) then
            if(bra_shell_1 == bra_shell_2.or.ket_shell_2 == bra_shell_2) then
                if(bra_shell_1 == ket_shell_2) call report_invalid_shell_transfer()
                if(ket_shell_2 /= bra_shell_2) then
                    call accumulate_single_transfer_triple_on_second_shell(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_1, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                else
                    call accumulate_single_transfer_triple_on_first_shell(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_1, &
                        & bra_shell_1,bra_shell_2,ket_shell_1, &
                        ket_shell_2,interaction_type)
                end if
            else
                call accumulate_single_transfer_with_spectator_shell(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_1, &
                    & bra_shell_2,2,bra_shell_1,bra_shell_2, &
                    ket_shell_2,ket_shell_1,interaction_type)
            end if
            return
        end if
        call report_invalid_shell_transfer()

    contains
        subroutine report_invalid_shell_transfer()
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN accumulate_two_electron_single_transfer  PMGG RAGG")')
            stop
        end subroutine report_invalid_shell_transfer
    end subroutine accumulate_two_electron_single_transfer
end module two_electron_single_transfer
