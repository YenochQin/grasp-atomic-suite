! Pack bra and ket quantum numbers for up to four active subshells.
module jj_active_subshell_states
    implicit none (type, external)

    ! Module interface
    private
    public :: perko2
contains
    subroutine perko2(first_shell_position,second_shell_position,third_shell_position,fourth_shell_position, &
        active_shell_count)
        ! Module dependencies
        use jj_tensor_coupling_tables, only: ket_shell_values_1 => bd1, ket_shell_values_2 => bd2, ket_shell_values_3 => bd3, &
            ket_shell_values_4 => bd4, bra_shell_values_1 => bk1, bra_shell_values_2 => bk2, &
            bra_shell_values_3 => bk3, bra_shell_values_4 => bk4, ket_shell_codes_1 => id1, ket_shell_codes_2 => id2, &
            ket_shell_codes_3 => id3, ket_shell_codes_4 => id4, bra_shell_codes_1 => ik1, bra_shell_codes_2 => ik2, &
            bra_shell_codes_3 => ik3, bra_shell_codes_4 => ik4
        use jj_subshell_state_packing, only: perko1
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_shell_position, second_shell_position, third_shell_position, &
            fourth_shell_position, active_shell_count
        call perko1(first_shell_position,bra_shell_values_1,bra_shell_codes_1,ket_shell_values_1,ket_shell_codes_1)
        if(active_shell_count == 1)return
        call perko1(second_shell_position,bra_shell_values_2,bra_shell_codes_2,ket_shell_values_2,ket_shell_codes_2)
        if(active_shell_count == 2)return
        call perko1(third_shell_position,bra_shell_values_3,bra_shell_codes_3,ket_shell_values_3,ket_shell_codes_3)
        if(active_shell_count == 3)return
        call perko1(fourth_shell_position,bra_shell_values_4,bra_shell_codes_4,ket_shell_values_4,ket_shell_codes_4)
        return
    end subroutine perko2
end module jj_active_subshell_states
