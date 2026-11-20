! Evaluate one creation or annihilation tensor in the few-electron high-j basis.
module jj_a_matrix_element_few_electrons
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_jj_a_few_electrons
contains
    subroutine evaluate_jj_a_few_electrons(bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values, &
        & operator_q_projection_1, &
        a_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, angular_tolerance
        implicit none (type, external)

        ! Arguments
        integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
        real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
        real(real64), intent(in)               :: operator_q_projection_1
        real(real64), intent(out)              :: a_matrix_element

        ! Local declarations
        integer      :: twice_phase_exponent, integer_norm_squared
        real(real64) :: norm_squared
        a_matrix_element=zero
        if(operator_q_projection_1 < angular_tolerance) then
            integer_norm_squared=(ket_shell_codes(6)+1)*ket_shell_codes(4)
            norm_squared=dble(integer_norm_squared)
            a_matrix_element=dsqrt(norm_squared)
            twice_phase_exponent=ket_shell_codes(6)+ket_shell_codes(3)-bra_shell_codes(6)+ket_shell_codes(4)*2
            if((twice_phase_exponent/4)*4 /= twice_phase_exponent)a_matrix_element=-a_matrix_element
        else
            integer_norm_squared=(bra_shell_codes(6)+1)*bra_shell_codes(4)
            norm_squared=dble(integer_norm_squared)
            a_matrix_element=dsqrt(norm_squared)
            if((bra_shell_codes(4)/2)*2 /= bra_shell_codes(4))a_matrix_element=-a_matrix_element
        end if
        return
    end subroutine evaluate_jj_a_few_electrons
end module jj_a_matrix_element_few_electrons
