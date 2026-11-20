module breit_radial_integral_1_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Charlotte Froese Fischer
        !                       Gediminas Gaigalas  10/05/17
        subroutine breit_radial_integral_1 (orbital_a,orbital_b,orbital_c,orbital_d,k,radial_integral_value)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_a, orbital_b, orbital_c, orbital_d, k
            real(real64), intent(out) :: radial_integral_value
        end subroutine breit_radial_integral_1
    end interface
end module breit_radial_integral_1_interface
