! Sort three distinct subshell positions for recoupling.
module three_subshell_ordering
    implicit none (type, external)

    ! Module interface
    private
    public :: eile
contains
    subroutine eile(shell_position_a,shell_position_b,shell_position_c,first_sorted_shell,middle_sorted_shell, &
        last_sorted_shell)
        implicit none (type, external)

        ! Arguments
        integer, intent(in)  :: shell_position_a, shell_position_b, shell_position_c
        integer, intent(out) :: first_sorted_shell, middle_sorted_shell, last_sorted_shell
        first_sorted_shell=shell_position_a
        last_sorted_shell=shell_position_a
        if(first_sorted_shell > shell_position_b)first_sorted_shell=shell_position_b
        if(last_sorted_shell < shell_position_b)last_sorted_shell=shell_position_b
        if(first_sorted_shell > shell_position_c)first_sorted_shell=shell_position_c
        if(last_sorted_shell < shell_position_c)last_sorted_shell=shell_position_c
        if((shell_position_a > first_sorted_shell).and.(shell_position_a < last_sorted_shell))middle_sorted_shell= &
            shell_position_a
        if((shell_position_b > first_sorted_shell).and.(shell_position_b < last_sorted_shell))middle_sorted_shell= &
            shell_position_b
        if((shell_position_c > first_sorted_shell).and.(shell_position_c < last_sorted_shell))middle_sorted_shell= &
            shell_position_c
        return
    end subroutine eile
end module three_subshell_ordering
