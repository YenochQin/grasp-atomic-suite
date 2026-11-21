! Dispatch a pair transfer that changes four distinct occupations.
module two_electron_four_shell_transfer
    implicit none (type, external)

    ! Module interface
    private
    public :: accumulate_four_shell_transfers
contains
    subroutine accumulate_four_shell_transfers(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
        & interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use csf_pair_angular_state,          only: peel_orbital_count
        use four_shell_transfer_12_34, only: accumulate_four_shell_transfer_12_34
        use four_shell_transfer_13_24, only: accumulate_four_shell_transfer_13_24
        use four_shell_transfer_14_23, only: accumulate_four_shell_transfer_14_23
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
            interaction_type
        if(peel_orbital_count <= 3)return
        ! Sort the four peel positions for the recoupling kernels. Case 1 means the first
        ! position belongs to the bra pair, and case 2 means it belongs to the ket pair.
        if(bra_shell_2 < ket_shell_1) then
            call accumulate_four_shell_transfer_12_34(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1, &
                & ket_shell_2,1,interaction_type)
        else if(bra_shell_1 > ket_shell_2.and.bra_shell_2 > ket_shell_2) then
            call accumulate_four_shell_transfer_12_34(bra_csf_index,ket_csf_index,ket_shell_1,ket_shell_2,bra_shell_1, &
                & bra_shell_2,2,interaction_type)
        else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 < ket_shell_2.and.bra_shell_1 < ket_shell_1) then
            call accumulate_four_shell_transfer_13_24(bra_csf_index,ket_csf_index,bra_shell_1,ket_shell_1,bra_shell_2, &
                & ket_shell_2,1,interaction_type)
        else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 > ket_shell_2.and.bra_shell_1 > ket_shell_1) then
            call accumulate_four_shell_transfer_13_24(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,ket_shell_2, &
                & bra_shell_2,2,interaction_type)
        else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 > ket_shell_2.and.bra_shell_1 < ket_shell_1) then
            call accumulate_four_shell_transfer_14_23(bra_csf_index,ket_csf_index,bra_shell_1,ket_shell_1,ket_shell_2, &
                & bra_shell_2,1,interaction_type)
        else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 < ket_shell_2.and.bra_shell_1 > ket_shell_1) then
            call accumulate_four_shell_transfer_14_23(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,bra_shell_2, &
                & ket_shell_2,2,interaction_type)
        else
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN accumulate_four_shell_transfers ")')
            stop
        end if
        return
    end subroutine accumulate_four_shell_transfers
end module two_electron_four_shell_transfer
