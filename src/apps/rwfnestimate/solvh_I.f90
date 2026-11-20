      module solvh_i
      implicit none (type, external)
      interface
      subroutine solvh (orbital_index, solve_failed)
      use grasp_parameters, only: max_radial_points, radial_workspace_capacity, max_orbitals
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      logical, intent(out) :: solve_failed
      end subroutine solvh
      end interface
      end module solvh_i
