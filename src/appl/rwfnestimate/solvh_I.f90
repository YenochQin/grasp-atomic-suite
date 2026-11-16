      module solvh_i
      implicit none (type, external)
      interface
      subroutine solvh (iorb, fail)
      use parameter_def, only: max_radial_points, radial_workspace_capacity, max_orbitals
      implicit none (type, external)
      integer, intent(in) :: iorb
      logical, intent(out) :: fail
      end subroutine
      end interface
      end module
