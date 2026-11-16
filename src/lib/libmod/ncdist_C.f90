      module ncdist_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,   only:  max_radial_points
      implicit none (type, external)
      private :: real64
      real(real64) :: zdist(max_radial_points)
      end module ncdist_c
