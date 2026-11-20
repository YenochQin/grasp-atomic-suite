module lookup_breit_radial_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        real(real64) function lookup_breit_radial_integral (itype, orbital_a, orbital_b, orbital_c, orbital_d, k)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: itype
            integer, intent(in) :: orbital_a
            integer, intent(in) :: orbital_b
            integer, intent(in) :: orbital_c
            integer, intent(in) :: orbital_d
            integer, intent(in) :: k
        end function lookup_breit_radial_integral
    end interface
end module lookup_breit_radial_integral_interface
