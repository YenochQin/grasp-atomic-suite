      module bessel_function_workspace
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only:  max_radial_points
      implicit none (type, external)
      private :: real64
      real(real64) :: wij(2)
      real(real64) :: bessj(2, 2, max_radial_points), bessn(2, 2, max_radial_points)
      real(real64) :: bj(max_radial_points, 3)
      real(real64) :: tc(max_radial_points), td(max_radial_points)
      end module bessel_function_workspace
