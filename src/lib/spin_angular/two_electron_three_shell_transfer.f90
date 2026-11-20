! Dispatch pair splitting or joining among three subshells.
module two_electron_three_shell_transfer
    implicit none (type, external)

    ! Module interface
    private
    public :: accumulate_three_shell_transfers
contains
    subroutine accumulate_three_shell_transfers(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
        & interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use csf_pair_angular_state,            only: peel_orbital_count
        use three_shell_pair_split_or_join, only: accumulate_three_shell_pair_split_or_join
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
            interaction_type
        if(peel_orbital_count <= 2)return
        if(bra_shell_1 == bra_shell_2) then
            call accumulate_three_shell_pair_split_or_join(bra_csf_index,ket_csf_index,ket_shell_1,ket_shell_2,bra_shell_1,1, &
                & bra_shell_1,bra_shell_2, &
                ket_shell_1,ket_shell_2,interaction_type)
        else if(ket_shell_1 == ket_shell_2) then
            call accumulate_three_shell_pair_split_or_join(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,2, &
                & bra_shell_1,bra_shell_2, &
                ket_shell_1,ket_shell_2,interaction_type)
        else
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN EL4 ")')
            stop
        end if
        return
    end subroutine accumulate_three_shell_transfers
end module two_electron_three_shell_transfer
