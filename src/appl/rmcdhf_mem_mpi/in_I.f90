      module in_i
      implicit none (type, external)
      interface
      subroutine in (iorb, jp, p, q, mtp)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_radial_points
      implicit none (type, external)
      integer, intent(in) :: iorb
      integer, intent(in) :: jp
      real(real64), intent(inout) :: p(max_radial_points)
      real(real64), intent(inout) :: q(max_radial_points)
      integer, intent(out) :: mtp
      end subroutine in
      end interface
      end module in_i
