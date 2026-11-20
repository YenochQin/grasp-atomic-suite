module energy_outside_bounds_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        logical function energy_outside_bounds (etry)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: etry
        end function energy_outside_bounds
    end interface
end module energy_outside_bounds_interface
