      module tail_i
      implicit none (type, external)
      interface
      subroutine tail (iorb, p, q, jp, mtp)
      use parameter_def, only: max_radial_points, radial_workspace_capacity, max_orbitals
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: iorb
      real(real64), intent(inout) :: p(max_radial_points)
      real(real64), intent(out) :: q(max_radial_points)
      integer, intent(in) :: jp
      integer, intent(out) :: mtp
      end subroutine
      end interface
      end module
