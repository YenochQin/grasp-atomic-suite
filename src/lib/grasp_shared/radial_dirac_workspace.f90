      module radial_dirac_workspace
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      implicit none (type, external)
      private :: real64
      integer :: mtp0
      real(real64) :: p(max_radial_points), q(max_radial_points)
      real(real64) :: p0, q0
      real(real64) :: tf(max_radial_points), tg(max_radial_points), xu(max_radial_points), xv(max_radial_points)
      end module radial_dirac_workspace
