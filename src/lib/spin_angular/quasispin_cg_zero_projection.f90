! Evaluate a rank-zero or rank-one quasispin coefficient with projection zero.
module quasispin_cg_zero_projection
    implicit none (type, external)

    ! Module interface
    private
    public :: cle0sm
contains
    subroutine cle0sm(ket_q,ket_q_projection,quasispin_rank,bra_q,bra_q_projection,cg_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, one, two, angular_tolerance => eps
        use angular_triangle_and_parity, only: ittk
        use quasispin_cg_rank_one_zero_projection, only: c1e0sm
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in)  :: ket_q, ket_q_projection, quasispin_rank, bra_q, bra_q_projection
        real(real64), intent(out) :: cg_coefficient

        ! Local declarations
        integer :: two_ket_q, two_bra_q, two_quasispin_rank
        cg_coefficient=zero
        two_ket_q=two*ket_q+tenth
        two_bra_q=two*bra_q+tenth
        two_quasispin_rank=two*quasispin_rank+tenth
        if(ittk(two_ket_q,two_bra_q,two_quasispin_rank)==0)return
        if(quasispin_rank<angular_tolerance) then
            if((ket_q+tenth)<dabs(ket_q_projection))return
            if((bra_q+tenth)<dabs(bra_q_projection))return
            if(dabs(ket_q-bra_q)>angular_tolerance)return
            if(dabs(ket_q_projection-bra_q_projection)>angular_tolerance)return
            cg_coefficient=one
        else
            call c1e0sm(ket_q,ket_q_projection,bra_q,bra_q_projection,cg_coefficient)
        end if
        return
    end subroutine cle0sm
end module quasispin_cg_zero_projection
