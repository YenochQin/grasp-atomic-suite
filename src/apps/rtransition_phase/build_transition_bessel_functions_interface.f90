module build_transition_bessel_functions_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_transition_bessel_functions (w)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: w
        end subroutine build_transition_bessel_functions
    end interface
end module build_transition_bessel_functions_interface
