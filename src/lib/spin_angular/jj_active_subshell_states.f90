! Pack bra and ket quantum numbers for up to four active subshells.
module jj_active_subshell_states
    implicit none (type, external)

    ! Module interface
    private
    public :: prepare_jj_active_subshell_states
contains
    subroutine prepare_jj_active_subshell_states(first_shell_position,second_shell_position,third_shell_position, &
        & fourth_shell_position, &
        active_shell_count)
        ! Module dependencies
        use jj_tensor_coupling_tables, only: ket_shell_values_1, ket_shell_values_2, ket_shell_values_3, ket_shell_values_4, &
            bra_shell_values_1, bra_shell_values_2, bra_shell_values_3, bra_shell_values_4, ket_shell_codes_1, &
            ket_shell_codes_2, ket_shell_codes_3, ket_shell_codes_4, bra_shell_codes_1, bra_shell_codes_2, bra_shell_codes_3, &
            bra_shell_codes_4
        use jj_subshell_state_packing, only: pack_jj_subshell_state
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_shell_position, second_shell_position, third_shell_position, &
            fourth_shell_position, active_shell_count
        call pack_jj_subshell_state(first_shell_position,bra_shell_values_1,bra_shell_codes_1,ket_shell_values_1,ket_shell_codes_1)
        if(active_shell_count == 1)return
        call pack_jj_subshell_state(second_shell_position,bra_shell_values_2,bra_shell_codes_2,ket_shell_values_2,ket_shell_codes_2)
        if(active_shell_count == 2)return
        call pack_jj_subshell_state(third_shell_position,bra_shell_values_3,bra_shell_codes_3,ket_shell_values_3,ket_shell_codes_3)
        if(active_shell_count == 3)return
        call pack_jj_subshell_state(fourth_shell_position,bra_shell_values_4,bra_shell_codes_4,ket_shell_values_4,ket_shell_codes_4)
        return
    end subroutine prepare_jj_active_subshell_states
end module jj_active_subshell_states
