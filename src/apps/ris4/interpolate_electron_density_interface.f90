module interpolate_electron_density_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine interpolate_electron_density (xa, ya, dens)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: xa(3)
            real(real64), intent(in) :: ya(3)
            real(real64), intent(out) :: dens
        end subroutine interpolate_electron_density
    end interface
end module interpolate_electron_density_interface
