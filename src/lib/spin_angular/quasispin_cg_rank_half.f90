! Evaluate a quasispin Clebsch–Gordan coefficient of rank one half.
module quasispin_cg_rank_half
    implicit none (type, external)
    private
    public :: c0t5s
contains
      subroutine c0t5s(ket_q,ket_q_projection,operator_q_projection,bra_q,bra_q_projection,cg_coefficient)
      use, intrinsic :: iso_fortran_env, only: real64
      use numeric_constants,          only: zero, tenth , half, one, two, angular_tolerance => eps
      use angular_triangle_and_parity, only: ittk
      implicit none (type, external)
      real(real64), intent(in)  :: ket_q, ket_q_projection, operator_q_projection, bra_q, bra_q_projection
      real(real64), intent(out) :: cg_coefficient
      integer                    :: two_ket_q, two_bra_q, projection_sign_index
      real(real64) :: projection_signs(2)
      projection_signs(1)=one
      projection_signs(2)=-one
      cg_coefficient=zero
      two_ket_q=two*ket_q+tenth
      two_bra_q=two*bra_q+tenth
      if(ittk(two_ket_q,two_bra_q,1) == 0)return
      if(dabs(ket_q_projection+operator_q_projection-bra_q_projection) > angular_tolerance)return
      if((half+tenth) < dabs(operator_q_projection))return
      if((ket_q+tenth) < dabs(ket_q_projection))return
      if((bra_q+tenth) < dabs(bra_q_projection))return
      projection_sign_index=dabs(half-operator_q_projection)+one+tenth
      if(dabs(ket_q+half-bra_q) < angular_tolerance) then
        cg_coefficient=dsqrt((bra_q+projection_signs(projection_sign_index)*bra_q_projection)/(two*bra_q))
      else
        if(dabs(ket_q-half-bra_q) > angular_tolerance)return
        cg_coefficient=-projection_signs(projection_sign_index)*dsqrt((bra_q- &
            projection_signs(projection_sign_index)*bra_q_projection+one)/(two*bra_q+two))
      end if
      return
      end subroutine c0t5s
end module quasispin_cg_rank_half
