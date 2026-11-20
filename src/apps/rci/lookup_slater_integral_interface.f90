module lookup_slater_integral_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine lookup_slater_integral (orbital_a, orbital_b, orbital_c, orbital_d, multipole_rank, tegral)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            !      INTEGER KMAX
            !      PARAMETER (KMAX = 20)
            integer, intent(inout) :: orbital_a
            integer, intent(inout) :: orbital_b
            integer, intent(inout) :: orbital_c
            integer, intent(inout) :: orbital_d
            integer, intent(in) :: multipole_rank
            real(real64), intent(out) :: tegral
        end subroutine lookup_slater_integral
    end interface
end module lookup_slater_integral_interface
