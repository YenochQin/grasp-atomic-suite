! Dispatch pair splitting or joining among three subshells.
module two_electron_three_shell_transfer
    implicit none (type, external)

    ! Module interface
    private
    public :: el4
contains
    subroutine el4(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,interaction_type)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use csf_pair_angular_state,            only: peel_orbital_count => npeel
        use three_shell_pair_split_or_join, only: el41
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
            interaction_type
        if(peel_orbital_count <= 2)return
        if(bra_shell_1 == bra_shell_2) then
            call el41(bra_csf_index,ket_csf_index,ket_shell_1,ket_shell_2,bra_shell_1,1,bra_shell_1,bra_shell_2, &
                ket_shell_1,ket_shell_2,interaction_type)
        else if(ket_shell_1 == ket_shell_2) then
            call el41(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,2,bra_shell_1,bra_shell_2, &
                ket_shell_1,ket_shell_2,interaction_type)
        else
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN EL4 ")')
            stop
        end if
        return
    end subroutine el4
end module two_electron_three_shell_transfer
