! Compute the angle of a Cartesian pair and map negative angles to the nonnegative branch.
module nonnegative_atan2
    implicit none (type, external)

    ! Module interface
    private
    public :: nonnegative_angle
contains
    real(real64) function nonnegative_angle (y_coordinate, x_coordinate)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: pi
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: y_coordinate
        real(real64), intent(in) :: x_coordinate

        ! Local declarations
        logical, save :: check_intrinsic_range = .true., intrinsic_returns_nonnegative_angle = .true.
        if (check_intrinsic_range) then
            nonnegative_angle = atan2(-1.0_real64,-1.0_real64)
            if (nonnegative_angle > 0.0_real64) then
                intrinsic_returns_nonnegative_angle = .true.
            else
                intrinsic_returns_nonnegative_angle = .false.
            end if
            check_intrinsic_range = .false.
        end if
        if (intrinsic_returns_nonnegative_angle) then
            nonnegative_angle = atan2(y_coordinate,x_coordinate)
        else
            if (y_coordinate >= 0.0_real64) then
                nonnegative_angle = atan2(y_coordinate,x_coordinate)
            else
                nonnegative_angle = pi + pi + atan2(y_coordinate,x_coordinate)
            end if
        end if
        return
    end function nonnegative_angle
end module nonnegative_atan2
