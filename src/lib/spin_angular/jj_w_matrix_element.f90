! Restore the occupation dependence of a pair tensor with quasispin coefficients.
module jj_w_matrix_element
    implicit none (type, external)

    ! Module interface
    private
    public :: wj1
contains
    subroutine wj1(bra_shell_codes,bra_shell_values,ket_shell_codes,ket_shell_values,w_angular_rank, &
        q_projection_1,q_projection_2,w_matrix_element)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, one, two, angular_tolerance => eps
        use subshell_tensor_diagnostics, only: mes
        use jj_w_matrix_element_high_j, only: w1jjg
        use quasispin_cg_zero_projection, only: cle0sm
        use quasispin_cg_rank_one_unit_projection, only: c1e1sm
        use jj_reduced_w, only: rwjj
        implicit none (type, external)

        ! Arguments
        integer,      intent(in)               :: w_angular_rank
        integer,      intent(in) :: bra_shell_codes(7), ket_shell_codes(7)
        real(real64), intent(in)               :: q_projection_1, q_projection_2
        real(real64), intent(in) :: bra_shell_values(3), ket_shell_values(3)
        real(real64), intent(out)              :: w_matrix_element

        ! Local declarations
        integer      :: w_quasispin_rank, net_q_projection, occupation_change_2
        real(real64) :: matrix_element_factor, combined_q_projection, reduced_w_element, real_quasispin_rank
        w_matrix_element=zero
        if(ket_shell_codes(3) == 9) then
            if(max0(bra_shell_codes(4),ket_shell_codes(4)) < 3) then
                if(bra_shell_codes(1) < 300) call mes(56)
                if(ket_shell_codes(1) < 300) call mes(56)
                occupation_change_2=q_projection_2+q_projection_2+q_projection_2*angular_tolerance
                if((ket_shell_codes(4)+occupation_change_2) > 2) call mes(2)
                call w1jjg(w_angular_rank,q_projection_1,q_projection_2,bra_shell_codes,bra_shell_values,ket_shell_codes, &
                    ket_shell_values,w_matrix_element)
                return
            else
                print*, "ERROR in  WJ1"
                stop
            end if
        else if(ket_shell_codes(3) > 9) then
            occupation_change_2=q_projection_2+q_projection_2+q_projection_2*angular_tolerance
            if((ket_shell_codes(4)+occupation_change_2) > 2) call mes(2)
            call w1jjg(w_angular_rank,q_projection_1,q_projection_2,bra_shell_codes,bra_shell_values,ket_shell_codes, &
                ket_shell_values,w_matrix_element)
            return
        end if

        ! The sum of the two quasispin projections fixes the net occupation change.
        combined_q_projection=q_projection_1+q_projection_2
        if(dabs(combined_q_projection) >= angular_tolerance) then
            if(((w_angular_rank/2)*2) /= w_angular_rank)return
            net_q_projection=combined_q_projection+combined_q_projection*tenth
            if((bra_shell_codes(4)-ket_shell_codes(4)-2*net_q_projection) /= 0)return
            call c1e1sm(ket_shell_values(1),ket_shell_values(3),combined_q_projection,bra_shell_values(1), &
                bra_shell_values(3),matrix_element_factor)
            if(dabs(matrix_element_factor) < angular_tolerance)return
            call rwjj(bra_shell_codes(3),bra_shell_codes(1),ket_shell_codes(1),1,w_angular_rank,reduced_w_element)
            w_matrix_element=matrix_element_factor*reduced_w_element/dsqrt(two*bra_shell_values(1)+one)
        else if(bra_shell_codes(4) /= ket_shell_codes(4))then
            return
        else if(w_angular_rank == 0) then
            if(ket_shell_codes(1) /= bra_shell_codes(1))return
            if(q_projection_1 < angular_tolerance) then
                matrix_element_factor=dble(ket_shell_codes(3)+1-ket_shell_codes(4))
            else
                matrix_element_factor=-dble(ket_shell_codes(4))
            end if
            w_matrix_element=matrix_element_factor*dsqrt(dble(bra_shell_codes(6)+1)/dble(bra_shell_codes(3)+1))
        else
            w_quasispin_rank=1
            if(((w_angular_rank/2)*2) /= w_angular_rank)w_quasispin_rank=0
            real_quasispin_rank=dble(w_quasispin_rank)
            call cle0sm(ket_shell_values(1),ket_shell_values(3),real_quasispin_rank,bra_shell_values(1), &
                bra_shell_values(3),matrix_element_factor)
            if(dabs(matrix_element_factor) < angular_tolerance)return
            call rwjj(bra_shell_codes(3),bra_shell_codes(1),ket_shell_codes(1),w_quasispin_rank,w_angular_rank, &
                reduced_w_element)
            matrix_element_factor=matrix_element_factor*reduced_w_element
            w_matrix_element=matrix_element_factor/dsqrt(two*two*bra_shell_values(1)+two)
            if(q_projection_1 >= angular_tolerance)return
            if(((w_angular_rank/2)*2) /= w_angular_rank)w_matrix_element=-w_matrix_element
        end if
        return
    end subroutine wj1
end module jj_w_matrix_element
