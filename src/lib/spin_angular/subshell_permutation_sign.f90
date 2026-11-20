! Count inversions of four active subshell positions to obtain the permutation sign.
module subshell_permutation_sign
    implicit none (type, external)

    ! Module interface
    private
    public :: jfaze
contains
    integer function jfaze(first_shell_position,second_shell_position,third_shell_position, &
        fourth_shell_position) result(permutation_sign)
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_shell_position, second_shell_position, third_shell_position, fourth_shell_position
        permutation_sign=1
        if(first_shell_position > second_shell_position)permutation_sign=-permutation_sign
        if(first_shell_position > third_shell_position)permutation_sign=-permutation_sign
        if(first_shell_position > fourth_shell_position)permutation_sign=-permutation_sign
        if(second_shell_position > third_shell_position)permutation_sign=-permutation_sign
        if(second_shell_position > fourth_shell_position)permutation_sign=-permutation_sign
        if(third_shell_position > fourth_shell_position)permutation_sign=-permutation_sign
        return
    end function jfaze
end module subshell_permutation_sign
