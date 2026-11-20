module one_electron_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine one_electron_integral (first_orbital,second_orbital,radial_integral_value)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(inout) :: first_orbital, second_orbital
            real(real64), intent(out) :: radial_integral_value
        end subroutine one_electron_integral
    end interface
end module one_electron_integral_interface
