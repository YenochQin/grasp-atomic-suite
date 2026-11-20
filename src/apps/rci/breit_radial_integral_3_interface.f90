module breit_radial_integral_3_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine breit_radial_integral_3 (orbital_a,orbital_b,orbital_c,orbital_d,radial_integral_order,radial_integral_value)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, radial_integral_order
            real(real64), intent(out) :: radial_integral_value
        end subroutine breit_radial_integral_3
    end interface
end module breit_radial_integral_3_interface
