! Test the range and parity of the bra and ket quasispin projections.
module quasispin_projection_selection
    implicit none (type, external)

    ! Module interface
    private
    public :: izas1
contains
    integer function izas1(two_ket_q,ket_q_projection,two_bra_q,bra_q_projection) result(projections_allowed)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use numeric_constants,          only: tenth, two
        implicit none (type, external)

        ! Arguments
        integer,      intent(in) :: two_ket_q, two_bra_q
        real(real64), intent(in) :: ket_q_projection, bra_q_projection

        ! Local declarations
        integer :: two_abs_ket_projection, two_abs_bra_projection
        projections_allowed=0
        two_abs_ket_projection=two*dabs(ket_q_projection)+tenth
        if(two_abs_ket_projection > two_ket_q)return
        if(mod(two_ket_q+two_abs_ket_projection,2) /= 0)return
        two_abs_bra_projection=two*dabs(bra_q_projection)+tenth
        if(two_abs_bra_projection > two_bra_q)return
        if(mod(two_bra_q+two_abs_bra_projection,2) /= 0)return
        projections_allowed=1
        return
    end  function izas1
end module quasispin_projection_selection
