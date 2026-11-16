      module sbstep_i
      implicit none (type, external)
      interface
      subroutine sbstep (iorb, nstrt, nend, p, q)
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def, only: max_radial_points, radial_workspace_capacity, max_orbitals
      implicit none (type, external)

      integer, intent(in) :: iorb
      integer, intent(in) :: nstrt
      integer, intent(in) :: nend
      real(real64), intent(inout) :: p(max_radial_points)
      real(real64), intent(inout) :: q(max_radial_points)
      end subroutine
      end interface
      end module
