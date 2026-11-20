      module hovlap_i
      implicit none (type, external)
      interface
      real(real64) function hovlap (p, q, mtpo, np, kappa, z)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      implicit none (type, external)
      real(real64), intent(in) :: p(max_radial_points)
      real(real64), intent(in) :: q(max_radial_points)
      integer, intent(in) :: mtpo
      integer :: np
      integer :: kappa
      real(real64) :: z
      end function hovlap
      end interface
      end module hovlap_i
