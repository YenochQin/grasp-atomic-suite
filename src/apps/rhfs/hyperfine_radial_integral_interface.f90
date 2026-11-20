module hyperfine_radial_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function hyperfine_radial_integral (first_orbital, second_orbital, k)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: first_orbital
            integer, intent(in) :: second_orbital
            integer, intent(in) :: k
        end function hyperfine_radial_integral
    end interface
end module hyperfine_radial_integral_interface
