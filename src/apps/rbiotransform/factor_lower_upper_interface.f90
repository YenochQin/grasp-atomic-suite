module factor_lower_upper_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine factor_lower_upper (a, l, u, ndim)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ndim
            real(real64), intent(in) :: a(ndim, ndim)
            real(real64), intent(out) :: l(*)
            real(real64), intent(inout) :: u(*)
        end subroutine factor_lower_upper
    end interface
end module factor_lower_upper_interface
