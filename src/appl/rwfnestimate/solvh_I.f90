      module solvh_i
      implicit none (type, external)
      interface
      subroutine solvh (iorb, fail)
      use grasp_parameters, only: max_radial_points, radial_workspace_capacity, max_orbitals
      implicit none (type, external)
      integer, intent(in) :: iorb
      logical, intent(out) :: fail
      end subroutine solvh
      end interface
      end module solvh_i
