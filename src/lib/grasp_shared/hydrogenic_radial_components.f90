      module hydrogenic_radial_components
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      implicit none (type, external)
      private :: real64
      real(real64) :: ph(max_radial_points), qh(max_radial_points)
      end module hydrogenic_radial_components
