      module scf_radial_potentials
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      implicit none (type, external)
      private :: real64
      real(real64) :: yp(max_radial_points), xp(max_radial_points), xq(max_radial_points)
      end module scf_radial_potentials
