! Compute the angle of a Cartesian pair and map negative angles to the nonnegative branch.
module nonnegative_atan2
    implicit none (type, external)
    private
    public :: arctan
contains
      real(real64) function arctan (y_coordinate, x_coordinate)
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: pi
      implicit none (type, external)
      real(real64), intent(in) :: y_coordinate
      real(real64), intent(in) :: x_coordinate
      logical, save :: check_intrinsic_range = .true., intrinsic_returns_nonnegative_angle = .true.
      if (check_intrinsic_range) then
         arctan = atan2(-1.0_real64,-1.0_real64)
         if (arctan > 0.0_real64) then
            intrinsic_returns_nonnegative_angle = .true.
         else
            intrinsic_returns_nonnegative_angle = .false.
         end if
         check_intrinsic_range = .false.
      end if
      if (intrinsic_returns_nonnegative_angle) then
         arctan = atan2(y_coordinate,x_coordinate)
      else
         if (y_coordinate >= 0.0_real64) then
            arctan = atan2(y_coordinate,x_coordinate)
         else
            arctan = pi + pi + atan2(y_coordinate,x_coordinate)
         end if
      end if
      return
      end function arctan
end module nonnegative_atan2
