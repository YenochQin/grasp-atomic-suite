module breit_radial_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function breit_radial_integral (itype, orbital_a, orbital_c, orbital_b, orbital_d, k)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: itype
            integer, intent(in) :: orbital_a
            integer, intent(in) :: orbital_c
            integer, intent(in) :: orbital_b
            integer, intent(in) :: orbital_d
            integer :: k
        end function breit_radial_integral
    end interface
end module breit_radial_integral_interface
