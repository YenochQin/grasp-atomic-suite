      module horb_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only: max_radial_points
      implicit none (type, external)
      private :: real64
      real(real64) :: ph(max_radial_points), qh(max_radial_points)
      end module horb_c
