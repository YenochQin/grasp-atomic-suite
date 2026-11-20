! Multiply two single-tensor subshell elements and their quasispin factors.
module subshell_product_a_a
    implicit none (type, external)

    ! Module interface
    private
    public :: gg12
contains
    subroutine gg12(bra_shell_codes_1,bra_shell_codes_2,bra_shell_values_1,bra_shell_values_2,ket_shell_codes_1, &
        ket_shell_codes_2,ket_shell_values_1,ket_shell_values_2,q_projection_1,q_projection_2, &
        subshell_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, angular_tolerance => eps
        use subshell_tensor_diagnostics, only: mes
        use quasispin_cg_rank_half, only: c0t5s
        use jj_reduced_cfp, only: rmeajj
        implicit none (type, external)

        ! Arguments
        integer,      intent(in) :: bra_shell_codes_1(7), bra_shell_codes_2(7), ket_shell_codes_1(7), &
            ket_shell_codes_2(7)
        real(real64), intent(in)               :: q_projection_1, q_projection_2
        real(real64), intent(in) :: bra_shell_values_1(3), bra_shell_values_2(3), ket_shell_values_1(3), &
            ket_shell_values_2(3)
        real(real64), intent(out)              :: subshell_matrix_element

        ! Local declarations
        integer :: occupation_change_1, occupation_change_2
        real(real64) :: cg_product, second_reduced_cfp, second_cg_factor, first_reduced_cfp
        subshell_matrix_element=zero
        if(bra_shell_codes_1(3) > 9) then
            if(bra_shell_codes_1(4) > 2) call mes(32)
            if(ket_shell_codes_1(4) > 2) call mes(32)
        end if
        if(bra_shell_codes_2(3) > 9) then
            if(bra_shell_codes_2(4) > 2) call mes(32)
            if(ket_shell_codes_2(4) > 2) call mes(32)
        end if
        occupation_change_1=q_projection_1+q_projection_1+tenth*q_projection_1
        if(bra_shell_codes_1(4) /= (ket_shell_codes_1(4)+occupation_change_1))return
        occupation_change_2=q_projection_2+q_projection_2+tenth*q_projection_2
        if(bra_shell_codes_2(4) /= (ket_shell_codes_2(4)+occupation_change_2))return
        call c0t5s(ket_shell_values_1(1),ket_shell_values_1(3),q_projection_1,bra_shell_values_1(1), &
            bra_shell_values_1(3),cg_product)
        if(dabs(cg_product) < angular_tolerance)return
        call c0t5s(ket_shell_values_2(1),ket_shell_values_2(3),q_projection_2,bra_shell_values_2(1), &
            bra_shell_values_2(3),second_cg_factor)
        if(dabs(second_cg_factor) < angular_tolerance)return
        cg_product=cg_product*second_cg_factor
        call rmeajj(bra_shell_codes_1(3),bra_shell_codes_1(1),bra_shell_codes_1(7),bra_shell_codes_1(6), &
            ket_shell_codes_1(1),ket_shell_codes_1(7),ket_shell_codes_1(6),first_reduced_cfp)
        if(dabs(first_reduced_cfp) < angular_tolerance)return
        call rmeajj(bra_shell_codes_2(3),bra_shell_codes_2(1),bra_shell_codes_2(7),bra_shell_codes_2(6), &
            ket_shell_codes_2(1),ket_shell_codes_2(7),ket_shell_codes_2(6),second_reduced_cfp)
        if(dabs(second_reduced_cfp) < angular_tolerance)return
        subshell_matrix_element=cg_product*first_reduced_cfp*second_reduced_cfp/dsqrt(dble((bra_shell_codes_1(7)+1)* &
            (bra_shell_codes_2(7)+1)))
        return
    end subroutine gg12
end module subshell_product_a_a
