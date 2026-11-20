module field_shift_radial_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function field_shift_radial_integral (first_orbital, second_orbital)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer :: first_orbital
            integer :: second_orbital
        end function field_shift_radial_integral
    end interface
end module field_shift_radial_integral_interface
