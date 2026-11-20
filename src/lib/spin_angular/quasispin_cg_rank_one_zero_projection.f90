! Evaluate a rank-one quasispin coefficient with projection zero.
module quasispin_cg_rank_one_zero_projection
    implicit none (type, external)

    ! Module interface
    private
    public :: c1e0sm
contains
    subroutine c1e0sm(ket_q,ket_q_projection,bra_q,bra_q_projection,cg_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, one, two, three, angular_tolerance => eps
        use angular_triangle_and_parity, only: ittk
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in)  :: ket_q, ket_q_projection, bra_q, bra_q_projection
        real(real64), intent(out) :: cg_coefficient

        ! Local declarations
        integer :: two_ket_q, two_bra_q, zero_projection_parity_sum, rank_one_coupling_case
        cg_coefficient=zero
        two_ket_q=two*ket_q+tenth
        two_bra_q=two*bra_q+tenth
        if(ittk(two_ket_q,two_bra_q,2) == 0)return
        if(dabs(ket_q_projection-bra_q_projection) > angular_tolerance) return
        if((ket_q+tenth) < dabs(ket_q_projection)) return
        if((bra_q+tenth) < dabs(bra_q_projection)) return
        if(dabs(ket_q_projection) <= angular_tolerance) then
            zero_projection_parity_sum=ket_q+bra_q+one+tenth
            if((zero_projection_parity_sum/2)*2 /= zero_projection_parity_sum) return
        end if
        rank_one_coupling_case=ket_q-bra_q+two+tenth
        if(rank_one_coupling_case <= 0) return
        if(rank_one_coupling_case > 3) return
        if (rank_one_coupling_case == 1) then
            cg_coefficient=dsqrt(((bra_q+bra_q_projection)*(bra_q-bra_q_projection))/((two*bra_q-one)*bra_q))
        else if (rank_one_coupling_case == 2) then
            cg_coefficient=bra_q_projection/dsqrt(bra_q*(bra_q+one))
        else if (rank_one_coupling_case == 3) then
            cg_coefficient=-dsqrt(((bra_q+bra_q_projection+one)*(bra_q-bra_q_projection+one))/((bra_q+one)*(two*bra_q+ &
                three)))
        end if
        return
    end subroutine c1e0sm
end module quasispin_cg_rank_one_zero_projection
