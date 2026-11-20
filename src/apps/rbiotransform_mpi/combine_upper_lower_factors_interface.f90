module combine_upper_lower_factors_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine combine_upper_lower_factors (a, u, l, ndim, scr)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ndim
            real(real64), intent(in) :: a(ndim, ndim)
            real(real64), intent(out) :: u(ndim, ndim)
            real(real64), intent(out) :: l(ndim, ndim)
            real(real64), intent(inout) :: scr(*)
        end subroutine combine_upper_lower_factors
    end interface
end module combine_upper_lower_factors_interface
