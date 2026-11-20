! Evaluate a rank-one quasispin coefficient with projection +1 or -1.
module quasispin_cg_rank_one_unit_projection
    implicit none (type, external)

    ! Module interface
    private
    public :: c1e1sm
contains
    subroutine c1e1sm(ket_q,ket_q_projection,operator_q_projection,bra_q,bra_q_projection,cg_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: zero, tenth, one, two, three, angular_tolerance => eps
        use angular_triangle_and_parity, only: ittk
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in)  :: ket_q, ket_q_projection, operator_q_projection, bra_q, bra_q_projection
        real(real64), intent(out) :: cg_coefficient

        ! Local declarations
        integer                    :: projection_sign_index, two_ket_q, two_bra_q
        real(real64) :: projection_signs(2)
        projection_signs(1)=one
        projection_signs(2)=-one
        cg_coefficient=zero
        two_ket_q=two*ket_q+tenth
        two_bra_q=two*bra_q+tenth
        if(ittk(two_ket_q,two_bra_q,2)==0)return
        if(dabs(ket_q_projection+operator_q_projection-bra_q_projection)>angular_tolerance)return
        if((ket_q+tenth)<dabs(ket_q_projection))return
        if((bra_q+tenth)<dabs(bra_q_projection))return
        projection_sign_index=0
        if(dabs(operator_q_projection-one)<angular_tolerance)projection_sign_index=1
        if(dabs(operator_q_projection+one)<angular_tolerance)projection_sign_index=2
        if(projection_sign_index==0)return
        if(dabs(ket_q+one-bra_q)<angular_tolerance) then
            cg_coefficient=dsqrt((bra_q+projection_signs(projection_sign_index)*bra_q_projection-one)*(bra_q+ &
                projection_signs(projection_sign_index)*bra_q_projection)/ &
                ((two*bra_q-one)*two*bra_q))
        else if(dabs(ket_q-bra_q)<angular_tolerance) then
            cg_coefficient=-projection_signs(projection_sign_index)*dsqrt((bra_q+ &
                projection_signs(projection_sign_index)*bra_q_projection)*(bra_q- &
                projection_signs(projection_sign_index)*bra_q_projection+one)/ &
                ((bra_q+one)*two*bra_q))
        else if(dabs(ket_q-one-bra_q)>angular_tolerance) then
            return
        else
            cg_coefficient=dsqrt((bra_q-projection_signs(projection_sign_index)*bra_q_projection+one)*(bra_q- &
                projection_signs(projection_sign_index)*bra_q_projection+two)/ &
                ((two*bra_q+two)*(two*bra_q+three)))
        end if
        return
    end subroutine c1e1sm
end module quasispin_cg_rank_one_unit_projection
