      module tail_i
      implicit none (type, external)
      interface
      subroutine tail (orbital_index, large_component_values, small_component_values, matching_grid_index, mtp)
      use grasp_parameters, only: max_radial_points, radial_workspace_capacity, max_orbitals
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      real(real64), intent(inout) :: large_component_values(max_radial_points)
      real(real64), intent(out) :: small_component_values(max_radial_points)
      integer, intent(in) :: matching_grid_index
      integer, intent(out) :: mtp
      end subroutine tail
      end interface
      end module tail_i
